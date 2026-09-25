"""Restart selection, immutable inputs, and output history boundaries."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
SPEC=importlib.util.spec_from_file_location('graphite_restart',ROOT/'slurry/drivers/gr_re2/run_graphite.py')
RUNNER=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(RUNNER)


def case():
    return json.loads((ROOT/'slurry/cases/pure_gr.json').read_text())


def metadata(config,step,rate=100):
    cfg,derived=RUNNER.resolve(config,rate)
    values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
    keys = RUNNER.SURFACE_CHECKPOINT_KEYS if cfg['interaction']['surface_adhesion'] else ('local_gap_fraction',)
    values['interaction_model_version'] = RUNNER.SURFACE_ADHESION_VERSION
    immutable = {key: values[key] for key in keys}
    immutable.update(shear_rate=rate, dt_s=derived['dt_s'], particle_count=cfg['particles']['count'], ranks=1)
    return {'engine':'pure_gr','format_version':1,'complete':True,'step':step,'ranks':1,
            'dt_s':derived['dt_s'],'shear_rate_s_inv':rate,
            'immutable_config':immutable,
            'files':{'state.bin':1,'initial_particles.csv':1,'lattice_rank_0.bin':1},
            'output_bytes':{'history.csv':9,'particles.csv':9}}


def fixture(root,step,rate=100,config=None):
    path=root/'checkpoints'/('checkpoint_%020d'%step);path.mkdir(parents=True)
    data=metadata(config or case(),step,rate)
    for name in data['files']:(path/name).write_bytes(b'x')
    RUNNER.write_json(path/'checkpoint.json',data)
    return path,data


class RestartTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory(prefix='restart-unit-')
        self.root=Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_periodic_defaults_and_invalid_retention(self):
        original=case();snapshot=copy.deepcopy(original)
        cfg,_=RUNNER.resolve(original)
        self.assertEqual(cfg['output']['checkpoint_every_steps'],0)
        self.assertEqual(cfg['output']['checkpoint_every_seconds'],350*60)
        self.assertEqual(cfg['output']['checkpoint_keep'],2)
        values=RUNNER.solver_values(cfg,self.root,self.root/'initial_particles.csv',8)
        self.assertEqual(values['checkpoint_every'],0)
        self.assertEqual(values['checkpoint_seconds'],350*60)
        self.assertEqual(original,snapshot)
        for key,value in [('checkpoint_every_steps',-1),('checkpoint_every_seconds',float('nan')),
                          ('checkpoint_keep',0)]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                broken=case();broken['output'][key]=value;RUNNER.resolve(broken)

    def test_latest_complete_ignores_partial_directory(self):
        fixture(self.root,4);new,_=fixture(self.root,8)
        partial=self.root/'checkpoints/checkpoint_00000000000000000012.partial'
        partial.mkdir();(partial/'checkpoint.json').write_text('{}')
        self.assertEqual(RUNNER.latest_checkpoint(self.root)[0],new)
        self.assertEqual(RUNNER.latest_checkpoint(new)[0],new)

    def test_old_run_has_clear_rejection(self):
        (self.root/'history.csv').write_text('step,time_s\n')
        with self.assertRaisesRegex(ValueError,'fluid state'):
            RUNNER.latest_checkpoint(self.root)

    def test_missing_rank_file_and_wrong_engine_are_rejected(self):
        path,data=fixture(self.root,8)
        (path/'lattice_rank_0.bin').unlink()
        with self.assertRaisesRegex(ValueError,'Missing or truncated'):
            RUNNER.read_checkpoint(path)
        data['engine']='gr_baseline';RUNNER.write_json(path/'checkpoint.json',data)
        with self.assertRaisesRegex(ValueError,'pure_gr'):
            RUNNER.read_checkpoint(path)

    def test_continuation_accepts_solver_and_output_changes(self):
        data=metadata(case(),8)
        changed=case();changed['numerics']['particle_max_iterations']+=10
        changed['output']['sample_every_steps']=1
        changed['output']['checkpoint_every_seconds']=60
        cfg,derived=RUNNER.resolve(changed,end_strain=20)
        self.assertTrue(RUNNER.validate_restart(cfg,derived,data,1))

    def test_changed_force_timestep_or_ranks_is_rejected(self):
        data=metadata(case(),8)
        for change in ('work','range','curvature_switch','curvature_cutoff','time','ranks'):
            with self.subTest(change=change),self.assertRaisesRegex(ValueError,'differ'):
                changed=case()
                if change=='work':changed['interaction']['adhesion_work_J_m2'] *= 0.5
                if change=='range':changed['interaction']['adhesion_range_m'] *= 0.5
                if change=='curvature_switch':changed['interaction']['curvature_switch_gap_m'] *= 1.1
                if change=='curvature_cutoff':changed['interaction']['curvature_cutoff_gap_m'] *= 1.1
                if change=='time':changed['numerics']['time_step_s']=1e-7
                cfg,derived=RUNNER.resolve(changed)
                RUNNER.validate_restart(cfg,derived,data,2 if change=='ranks' else 1)

    def test_model_migration_is_rejected_in_both_directions(self):
        legacy = case()
        legacy['interaction'].pop('surface_adhesion')
        for key in RUNNER.SURFACE_ADHESION_DEFAULTS:
            legacy['interaction'].pop(key)
        legacy['interaction'].update(RUNNER.LOCAL_ADHESION_DEFAULTS)
        legacy['interaction']['local_gap_fraction'] = 1.0
        for old, new in ((legacy, case()), (case(), legacy)):
            with self.subTest(old_surface=old['interaction'].get('surface_adhesion',False)):
                cfg, derived = RUNNER.resolve(new)
                with self.assertRaisesRegex(ValueError, 'interaction law differs'):
                    RUNNER.validate_restart(cfg, derived, metadata(old, 8), 1)
        cfg, derived = RUNNER.resolve(legacy)
        self.assertTrue(RUNNER.validate_restart(cfg, derived, metadata(legacy, 8), 1))

    def test_surface_fingerprint_cannot_silently_omit_parameters(self):
        cfg, derived = RUNNER.resolve(case())
        for key in RUNNER.SURFACE_CHECKPOINT_KEYS:
            saved = metadata(case(), 8)
            saved['immutable_config'].pop(key)
            with self.subTest(key=key), self.assertRaises(ValueError):
                RUNNER.validate_restart(cfg, derived, saved, 1)
        saved = metadata(case(), 8)
        saved['immutable_config']['interaction_model_version'] += 1
        with self.assertRaisesRegex(ValueError, 'version/fingerprint'):
            RUNNER.validate_restart(cfg, derived, saved, 1)

    def test_endpoint_is_cumulative_and_complete_can_be_skipped(self):
        cfg,derived=RUNNER.resolve(case());data=metadata(case(),16)
        with self.assertRaisesRegex(ValueError,'endpoint'):
            RUNNER.validate_restart(cfg,derived,data,1,16)
        self.assertFalse(RUNNER.validate_restart(cfg,derived,data,1,16,allow_complete=True))
        self.assertTrue(RUNNER.validate_restart(cfg,derived,data,1,17))

    def test_history_copies_only_committed_prefix(self):
        path,data=fixture(self.root,8)
        for name in ('history.csv','particles.csv'):(self.root/name).write_bytes(b'header\n8\n12\n')
        output=self.root/'new';output.mkdir()
        self.assertTrue(RUNNER.copy_history_prefix(path,data,output))
        for name in ('history.csv','particles.csv'):
            self.assertEqual((output/name).read_bytes(),b'header\n8\n')
            self.assertEqual((self.root/name).read_bytes(),b'header\n8\n12\n')

    def test_history_truncation_is_rejected_and_missing_history_is_optional(self):
        path,data=fixture(self.root,8);output=self.root/'new';output.mkdir()
        self.assertFalse(RUNNER.copy_history_prefix(path,data,output))
        for name in ('history.csv','particles.csv'):(self.root/name).write_bytes(b'x')
        with self.assertRaisesRegex(ValueError,'shorter'):
            RUNNER.copy_history_prefix(path,data,output)

    def test_batch_restart_retains_rates_and_skips_finished_case(self):
        config=case();config['flow']['shear_rate_s_inv']=333
        config_path=self.root/'current.json';RUNNER.write_json(config_path,config)
        fixture(self.root/'pure_gr/g000_100',16,100,config)
        pending,_=fixture(self.root/'pure_gr/g001_10',8,10,config)
        result=subprocess.run([sys.executable,str(ROOT/'slurry/tools/run_slurry.py'),
                               '--restart',str(self.root),'--config',str(config_path),'--ranks','1',
                               '--max-steps','16','--dry-run'],stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,universal_newlines=True,check=True)
        plan=json.loads(result.stdout)
        self.assertEqual(len(plan['jobs']),1)
        self.assertEqual(len(plan['restart_skipped_complete']),1)
        job=plan['jobs'][0];args=job['argv']
        self.assertEqual(job['engine'],'pure_gr')
        self.assertEqual(float(args[args.index('--shear-rate')+1]),10)
        self.assertEqual(args[args.index('--restart')+1],str(pending))
        self.assertEqual(args[args.index('--config')+1],str(config_path))


if __name__=='__main__':unittest.main()
