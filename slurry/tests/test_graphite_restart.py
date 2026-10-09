"""Restart selection, immutable inputs, and output history boundaries."""
import copy
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
SPEC=importlib.util.spec_from_file_location('graphite_restart',ROOT/'slurry/drivers/gr_re2/run_graphite.py')
RUNNER=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(RUNNER)


def case(name='pure_gr'):
    return json.loads((ROOT/'slurry/cases'/('gr_CMC.json' if name=='gr_cmc' else 'pure_gr.json')).read_text())


def metadata(config,step,rate=100):
    cfg,derived=RUNNER.resolve(config,rate)
    values = RUNNER.solver_values(cfg, Path('run'), Path('particles.csv'), 0)
    keys = RUNNER.SURFACE_CHECKPOINT_KEYS if cfg['interaction']['surface_adhesion'] else ('local_gap_fraction',)
    values['interaction_model_version'] = RUNNER.SURFACE_ADHESION_VERSION
    immutable = {key: values[key] for key in keys}
    for physical,version_key,version,checkpoint_keys in (
            ('free_cmc_repulsion_pressure','free_cmc_repulsion_version',RUNNER.FREE_CMC_REPULSION_VERSION,RUNNER.FREE_CMC_CHECKPOINT_KEYS),
            ('free_cmc_inner_repulsion_work','free_cmc_inner_repulsion_version',RUNNER.FREE_CMC_INNER_REPULSION_VERSION,RUNNER.INNER_CMC_CHECKPOINT_KEYS),
            ('cmc_net_blend','cmc_net_potential_version',RUNNER.CMC_NET_POTENTIAL_VERSION,RUNNER.NET_CMC_CHECKPOINT_KEYS),
            ('cmc_coordination_enabled','cmc_coordination_version',RUNNER.CMC_COORDINATION_VERSION,RUNNER.CMC_COORDINATION_CHECKPOINT_KEYS),
            ('current_adhesion_rolling','current_adhesion_rolling_version',RUNNER.CURRENT_ADHESION_ROLLING_VERSION,RUNNER.CURRENT_ADHESION_ROLLING_CHECKPOINT_KEYS)):
        if values.get(physical,0)>0:
            values[version_key]=version
            immutable.update({key:values[key] for key in checkpoint_keys})
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
        self.environment={key:value for key,value in os.environ.items() if not key.startswith('SLURM_')}

    def tearDown(self):
        self.temporary.cleanup()

    def controller(self,*args,**kwargs):
        return subprocess.run([sys.executable,str(ROOT/'slurry/tools/run_slurry.py')]+list(args),
                              cwd=str(ROOT),env=kwargs.get('env',self.environment),
                              stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)

    def snapshot(self,config,engine='gr_cmc'):
        path=self.root/'input/cases'/(engine+'.json');path.parent.mkdir(parents=True,exist_ok=True)
        RUNNER.write_json(path,config)
        return path

    def contents(self):
        return {str(path.relative_to(self.root)):path.read_bytes() for path in self.root.rglob('*') if path.is_file()}

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

    def test_restart_infers_saved_cmc_and_preserves_original_rate_directories(self):
        config=case('gr_cmc')
        config['cmc']['net_potential']['contact_force_N']*=.37
        saved=self.snapshot(config)
        fixture(self.root/'gr_cmc/g009_10',8,10,config)
        fixture(self.root/'gr_cmc/g003_100',16,100,config)
        before=self.contents()
        result=self.controller('--restart',str(self.root),'--max-steps','32','--dry-run')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        plan=json.loads(result.stdout)
        self.assertEqual(Path(plan['output']),self.root)
        self.assertEqual(plan['ranks'],{'gr_cmc':1})
        self.assertEqual(len(plan['jobs']),2)
        observed=set()
        for job in plan['jobs']:
            self.assertEqual(job['engine'],'gr_cmc')
            argv=job['argv']
            self.assertEqual(Path(argv[1]),ROOT/'slurry/drivers/gr_re2/run_graphite.py')
            self.assertEqual(Path(argv[argv.index('--config')+1]),saved)
            run=Path(argv[argv.index('--output')+1]);observed.add(run.name)
            checkpoint=Path(argv[argv.index('--restart')+1])
            self.assertEqual(checkpoint.parent.parent,run)
            resolved=subprocess.run(argv+['--dry-run'],cwd=str(ROOT),env=self.environment,
                                    stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
            self.assertEqual(resolved.returncode,0,resolved.stderr)
            actual=json.loads(resolved.stdout)['config']['cmc']['net_potential']['contact_force_N']
            self.assertEqual(actual,config['cmc']['net_potential']['contact_force_N'])
        self.assertEqual(observed,{'g009_10','g003_100'})
        self.assertEqual(self.contents(),before,'dry-run must not touch saved output or inputs')

    def test_saved_input_is_authoritative_and_endpoint_override_is_cumulative(self):
        config=case('gr_cmc');cfg,derived=RUNNER.resolve(config,10,end_strain=1)
        step=derived['steps_to_requested_end_strain']
        directory=self.root/'gr_cmc/g000_10'
        fixture(directory,step,10,config)
        # The old effective configuration deliberately differs. Restart should
        # read the input snapshot that the user edited to extend the run.
        stale=copy.deepcopy(cfg);stale['cmc']['net_potential']['contact_force_N']*=2
        RUNNER.write_json(directory/'effective_config.json',stale)
        config['flow']['end_strain']=1
        saved=self.snapshot(config)
        result=self.controller('--restart',str(self.root),'--end-strain','2','--dry-run')
        self.assertEqual(result.returncode,0,result.stderr)
        argv=json.loads(result.stdout)['jobs'][0]['argv']
        self.assertEqual(Path(argv[argv.index('--config')+1]),saved)
        self.assertEqual(float(argv[argv.index('--end-strain')+1]),2.)
        self.assertEqual(json.loads(saved.read_text())['flow']['end_strain'],1)

    def test_restart_uses_checkpoint_ranks_in_larger_new_allocation(self):
        config=case('gr_cmc');self.snapshot(config)
        fixture(self.root/'gr_cmc/g000_100',8,100,config)
        result=self.controller('--restart',str(self.root),'--max-steps','16','--dry-run',
                               env=dict(self.environment,SLURM_NTASKS='8'))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(json.loads(result.stdout)['ranks'],{'gr_cmc':1})
        invalid=self.controller('--restart',str(self.root),'--max-steps','16','--ranks','2','--dry-run')
        self.assertNotEqual(invalid.returncode,0)
        self.assertIn('differ',invalid.stderr)

    def test_missing_snapshot_uses_rate_effective_config_not_repository_default(self):
        config=case('gr_cmc');config['cmc']['net_potential']['contact_force_N']*=.4
        directory=self.root/'gr_cmc/g007_100';fixture(directory,8,100,config)
        saved=directory/'effective_config.json';RUNNER.write_json(saved,RUNNER.resolve(config,100)[0])
        result=self.controller('--restart',str(self.root),'--max-steps','16','--dry-run')
        self.assertEqual(result.returncode,0,result.stderr)
        argv=json.loads(result.stdout)['jobs'][0]['argv']
        self.assertEqual(Path(argv[argv.index('--config')+1]),saved)
        saved.unlink()
        before=self.contents()
        missing=self.controller('--restart',str(self.root),'--max-steps','16','--dry-run')
        self.assertNotEqual(missing.returncode,0,'missing saved inputs must not silently select current repository physics')
        self.assertEqual(self.contents(),before)

    def test_explicit_changed_cmc_physics_rejected_without_modifying_run(self):
        config=case('gr_cmc');self.snapshot(config)
        fixture(self.root/'gr_cmc/g000_100',8,100,config)
        changed=copy.deepcopy(config);changed['cmc']['net_potential']['contact_force_N']*=2
        override=self.root/'changed.json';RUNNER.write_json(override,changed)
        before=self.contents()
        result=self.controller('--restart',str(self.root),'--config',str(override),
                               '--max-steps','16','--dry-run')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('differ',result.stderr)
        self.assertEqual(self.contents(),before)

    def test_saved_pure_gr_restart_remains_pure_gr(self):
        config=case();saved=self.snapshot(config,'pure_gr')
        directory=self.root/'pure_gr/g004_10';checkpoint,_=fixture(directory,8,10,config)
        result=self.controller('--restart',str(self.root),'--max-steps','16','--dry-run')
        self.assertEqual(result.returncode,0,result.stderr)
        job=json.loads(result.stdout)['jobs'][0];argv=job['argv']
        self.assertEqual(job['engine'],'pure_gr')
        self.assertEqual(Path(argv[argv.index('--config')+1]),saved)
        self.assertEqual(Path(argv[argv.index('--restart')+1]),checkpoint)
        self.assertEqual(Path(argv[argv.index('--output')+1]),directory)

    def test_restart_does_not_read_current_repository_run_settings(self):
        config=case('gr_cmc');self.snapshot(config)
        fixture(self.root/'gr_cmc/g000_100',8,100,config)
        result=self.controller('--restart',str(self.root),'--max-steps','16','--dry-run',
                               '--settings',str(self.root/'missing_current_settings.json'))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(json.loads(result.stdout)['jobs'][0]['engine'],'gr_cmc')

    def test_restart_skips_saved_case_that_never_wrote_a_checkpoint(self):
        self.snapshot(case(),'pure_gr');config=case('gr_cmc');self.snapshot(config)
        (self.root/'pure_gr').mkdir()
        fixture(self.root/'gr_cmc/g000_100',8,100,config)
        result=self.controller('--restart',str(self.root),'--max-steps','16','--dry-run')
        self.assertEqual(result.returncode,0,result.stderr)
        plan=json.loads(result.stdout)
        self.assertEqual([job['engine'] for job in plan['jobs']],['gr_cmc'])
        self.assertEqual(plan['restart_skipped_unstarted'],[str(self.root/'pure_gr')])

    def test_archived_cli_numerics_restored_but_old_endpoint_not_replayed(self):
        config=case('gr_cmc');self.snapshot(config)
        actual=copy.deepcopy(config);actual['numerics']['target_mach']*=.5
        directory=self.root/'gr_cmc/g000_100';fixture(directory,8,100,actual)
        archive=self.root/'restart_attempts/old';archive.mkdir(parents=True)
        RUNNER.write_json(archive/'batch_manifest.json',{
            'output':str(self.root),'jobs':[{'engine':'gr_cmc','input_numerics':config['numerics'],
                'argv':['python','old_driver','--output',str(directory),
                        '--target-mach',str(actual['numerics']['target_mach']),
                        '--end-strain','0.0001','--max-steps','8']}]})
        RUNNER.write_json(self.root/'batch_manifest.json',{'output':str(self.root),'jobs':[]})
        result=self.controller('--restart',str(self.root),'--max-steps','16','--dry-run')
        self.assertEqual(result.returncode,0,result.stderr)
        argv=json.loads(result.stdout)['jobs'][0]['argv']
        self.assertEqual(float(argv[argv.index('--target-mach')+1]),actual['numerics']['target_mach'])
        self.assertEqual(argv[argv.index('--max-steps')+1],'16')
        self.assertNotIn('--end-strain',argv)

    def test_concurrent_controller_rejected_before_run_mutation(self):
        config=case('gr_cmc');self.snapshot(config)
        fixture(self.root/'gr_cmc/g000_100',8,100,config)
        with (self.root/'.run_slurry.lock').open('a') as lock:
            fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
            before=self.contents()
            result=self.controller('--restart',str(self.root),'--max-steps','16')
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Another batch controller',result.stderr)
            self.assertEqual(self.contents(),before)

    def test_explicit_case_cannot_relabel_another_case_rate_directory(self):
        config=case('gr_cmc');self.snapshot(config)
        directory=self.root/'gr_cmc/g000_100';fixture(directory,8,100,config)
        result=self.controller('--restart',str(directory),'--cases','pure_gr',
                               '--max-steps','16','--dry-run')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('gr_cmc',result.stderr)


if __name__=='__main__':unittest.main()
