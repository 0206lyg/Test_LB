"""Common controller routing and isolated restarts for CMC-coated graphite."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
CONTROLLER = ROOT / 'slurry/tools/run_slurry.py'
SPEC = importlib.util.spec_from_file_location(
    'gr_cmc_controller_driver', ROOT / 'slurry/drivers/gr_re2/run_graphite.py')
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)
COMMON_SPEC = importlib.util.spec_from_file_location(
    'gr_cmc_controller_common', ROOT / 'slurry/tools/common.py')
COMMON = importlib.util.module_from_spec(COMMON_SPEC)
COMMON_SPEC.loader.exec_module(COMMON)


def checkpoint(directory, config, step, rate):
    cfg, derived = DRIVER.resolve(config, rate)
    values = DRIVER.solver_values(cfg, directory, directory / 'particles.csv', 0)
    values['interaction_model_version'] = DRIVER.SURFACE_ADHESION_VERSION
    immutable = {key: values[key] for key in DRIVER.SURFACE_CHECKPOINT_KEYS}
    if values.get('free_cmc_repulsion_pressure', 0) > 0:
        values['free_cmc_repulsion_version'] = DRIVER.FREE_CMC_REPULSION_VERSION
        immutable.update({key: values[key] for key in DRIVER.FREE_CMC_CHECKPOINT_KEYS})
    if values.get('free_cmc_inner_repulsion_work', 0) > 0:
        values['free_cmc_inner_repulsion_version'] = DRIVER.FREE_CMC_INNER_REPULSION_VERSION
        immutable.update({key: values[key] for key in DRIVER.INNER_CMC_CHECKPOINT_KEYS})
    if values.get('cmc_net_blend', 0) > 0:
        values['cmc_net_potential_version'] = DRIVER.CMC_NET_POTENTIAL_VERSION
        immutable.update({key: values[key] for key in DRIVER.NET_CMC_CHECKPOINT_KEYS})
    if values.get('cmc_coordination_enabled', 0) > 0:
        values['cmc_coordination_version'] = DRIVER.CMC_COORDINATION_VERSION
        immutable.update({key: values[key] for key in DRIVER.CMC_COORDINATION_CHECKPOINT_KEYS})
    if values.get('current_adhesion_rolling', 0) > 0:
        values['current_adhesion_rolling_version'] = DRIVER.CURRENT_ADHESION_ROLLING_VERSION
        immutable.update({key: values[key] for key in DRIVER.CURRENT_ADHESION_ROLLING_CHECKPOINT_KEYS})
    immutable.update(shear_rate=rate, dt_s=derived['dt_s'],
                     particle_count=cfg['particles']['count'], ranks=1)
    saved = directory / 'checkpoints' / ('checkpoint_%020d' % step)
    saved.mkdir(parents=True)
    files = {'state.bin': 1, 'initial_particles.csv': 1, 'lattice_rank_0.bin': 1}
    for name in files:
        (saved / name).write_bytes(b'x')
    DRIVER.write_json(saved / 'checkpoint.json', {
        'engine': 'pure_gr', 'format_version': 1, 'complete': True,
        'step': step, 'ranks': 1, 'dt_s': derived['dt_s'],
        'shear_rate_s_inv': rate, 'immutable_config': immutable,
        'files': files, 'output_bytes': {'history.csv': 0, 'particles.csv': 0}})
    return saved


class GrCmcControllerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='gr-cmc-controller-')
        self.root = Path(self.temporary.name)
        self.environment = {key: value for key, value in os.environ.items()
                            if not key.startswith('SLURM_')}

    def tearDown(self):
        self.temporary.cleanup()

    def run_controller(self, *args, **kwargs):
        output=[] if '--restart' in args else ['--output',str(self.root/'output')]
        return subprocess.run(
            [sys.executable, str(CONTROLLER), '--dry-run'] + output + list(args),
            cwd=str(ROOT), env=kwargs.get('env', self.environment),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True)

    def plan(self, *args, **kwargs):
        result = self.run_controller(*args, **kwargs)
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout)

    def test_default_cases_remain_cmc_and_bare_graphite(self):
        plan = self.plan()
        self.assertEqual([job['engine'] for job in plan['jobs']], ['pure_cmc', 'pure_gr'])

    def test_gr_cmc_re2_options_reach_the_shared_driver(self):
        environment = dict(self.environment, SLURM_NTASKS='3')
        plan = self.plan('--cases', 'gr_cmc', '--shear-rates', '10,100',
                         '--max-steps', '16', '--target-mach', '0.04',
                         '--time-step', '0', '--end-strain', '20', env=environment)
        self.assertEqual(plan['ranks'], {'gr_cmc': 3})
        self.assertEqual(len(plan['jobs']), 2)
        for job in plan['jobs']:
            args = job['argv']
            self.assertEqual(job['engine'], 'gr_cmc')
            self.assertEqual(Path(args[1]), ROOT / 'slurry/drivers/gr_re2/run_graphite.py')
            self.assertEqual(Path(args[args.index('--config') + 1]),
                             ROOT / 'slurry/cases/gr_CMC.json')
            self.assertIn('/gr_cmc/', args[args.index('--output') + 1])
            for option, value in (('--max-steps', '16'), ('--target-mach', '0.04'),
                                  ('--time-step', '0.0'), ('--end-strain', '20.0')):
                self.assertEqual(args[args.index(option) + 1], value)
            self.assertNotIn('--benchmark-steps', args)
            result = subprocess.run(args + ['--dry-run'], cwd=str(ROOT),
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, env=environment)
            self.assertEqual(result.returncode, 0, result.stderr)
            resolved = json.loads(result.stdout)
            self.assertEqual(resolved['config']['flow']['end_strain'], 20.0)
            self.assertGreater(resolved['derived']['dt_s'], 0)

    def test_relative_config_override_and_smoke_use_re2_step_cap(self):
        plan = self.plan('--cases', 'gr_cmc', '--config', 'slurry/cases/gr_CMC.json', '--smoke')
        args = plan['jobs'][0]['argv']
        self.assertEqual(Path(args[args.index('--config') + 1]), ROOT / 'slurry/cases/gr_CMC.json')
        self.assertEqual(args[args.index('--max-steps') + 1], '20')
        self.assertEqual(float(args[args.index('--shear-rate') + 1]), 100)

    def test_baseline_still_rejects_re2_numerical_options(self):
        for option, value in (('--target-mach', '0.04'), ('--time-step', '0'), ('--end-strain', '20')):
            with self.subTest(option=option):
                result = self.run_controller('--cases', 'gr_baseline', option, value)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('require pure_gr or gr_cmc', result.stderr)

    def test_restart_filters_mixed_engine_folder_and_skips_complete_rates(self):
        config = json.loads((ROOT / 'slurry/cases/gr_CMC.json').read_text())
        saved=self.root/'input/cases/gr_cmc.json';saved.parent.mkdir(parents=True)
        COMMON.write_json(saved,config)
        checkpoint(self.root / 'pure_gr/g000_999', config, 8, 999)
        checkpoint(self.root / 'gr_cmc/g000_100', config, 16, 100)
        pending = checkpoint(self.root / 'gr_cmc/g001_10', config, 8, 10)
        plan = self.plan('--cases', 'gr_cmc', '--restart', str(self.root),
                         '--max-steps', '16', '--ranks', '1')
        self.assertEqual(len(plan['jobs']), 1)
        self.assertEqual(len(plan['restart_skipped_complete']), 1)
        job = plan['jobs'][0]
        args = job['argv']
        self.assertEqual(job['engine'], 'gr_cmc')
        self.assertEqual(args[args.index('--restart') + 1], str(pending))
        self.assertEqual(float(args[args.index('--shear-rate') + 1]), 10)

    def test_restart_does_not_fall_back_to_other_engine(self):
        config = json.loads((ROOT / 'slurry/cases/pure_gr.json').read_text())
        checkpoint(self.root / 'pure_gr/g000_100', config, 8, 100)
        result = self.run_controller('--cases', 'gr_cmc', '--restart', str(self.root), '--max-steps', '16')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('gr_cmc', result.stderr)

    def test_actual_controller_snapshots_and_resolves_cmc_for_solver(self):
        # Run the real controller, copied driver, and generator; replace only
        # the HPC executable and use one particle to keep preparation cheap.
        config = json.loads((ROOT / 'slurry/cases/gr_CMC.json').read_text())
        config['particles'].update(count=1, initialization_mc_sweeps=0)
        config_path = self.root / 'small_gr_cmc.json'
        COMMON.write_json(config_path, config)
        executable = self.root / 'fake_slurry'
        build_info = {'mpi_enabled': False, 'rough_contact': True,
                      'local_gap_adhesion': True, 'surface_adhesion_version': 1,
                      'pass_max_version': 1, 'cmc_contact_version': 1,
                      'free_cmc_repulsion_version': 1, 'free_cmc_inner_repulsion_version': 1,
                      'cmc_net_potential_version': 1,
                      'cmc_coordination_version': 1,
                      'current_adhesion_rolling_version': 1,
                      'pure_gr_checkpoint_version': 1, 'particle_solver': 'petsc',
                      'engines': ['pure_gr', 'gr_cmc']}
        executable.write_text(
            '#!' + sys.executable + '\n'
            'import json, os, pathlib, sys\n'
            'if sys.argv[1:] == ["--build-info"]:\n'
            '    print(' + repr(json.dumps(build_info)) + ')\n'
            'else:\n'
            '    pathlib.Path("observed_solver.json").write_text(json.dumps({\n'
            '        "engine": os.environ.get("SLURRY_ENGINE"), "argv": sys.argv[1:]}))\n'
            '    pathlib.Path("status.json").write_text(json.dumps({"status": "MAX_STEPS"}))\n')
        executable.chmod(0o755)
        COMMON.write_json(self.root / 'build_manifest.json', {
            'build_id': 'controller-fixture', 'build_info': build_info,
            'executable_sha256': COMMON.digest(executable),
            'source_inputs': COMMON.source_inputs()})
        output = self.root / 'actual_output'
        result = subprocess.run(
            [sys.executable, str(CONTROLLER), '--cases', 'gr_cmc',
             '--config', str(config_path), '--executable', str(executable),
             '--output', str(output), '--shear-rates', '100', '--max-steps', '2',
             '--end-strain', '20', '--target-mach', '0.04', '--ranks', '1'],
            cwd=str(ROOT), env=self.environment, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, universal_newlines=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(json.loads((output / 'batch_status.json').read_text())['status'], 'COMPLETED')
        batch = json.loads((output / 'batch_manifest.json').read_text())
        args = batch['jobs'][0]['argv']
        self.assertEqual(Path(args[1]), output / 'input/drivers/gr_re2/run_graphite.py')
        self.assertEqual(Path(args[args.index('--config') + 1]), output / 'input/cases/gr_cmc.json')
        self.assertEqual(json.loads((output / 'input/cases/gr_cmc.json').read_text()), config)
        run = output / 'gr_cmc/g000_100'
        observed = json.loads((run / 'observed_solver.json').read_text())
        self.assertEqual(observed['engine'], 'gr_cmc')
        self.assertEqual(observed['argv'], ['--config', str(run / 'resolved_run.cfg')])
        resolved = json.loads((run / 'effective_config.json').read_text())
        self.assertEqual(resolved['cmc']['q_sat'], config['cmc']['q_sat'])
        self.assertEqual(resolved['flow']['end_strain'], 20)
        self.assertEqual(resolved['numerics']['target_mach'], 0.04)
        values = dict(line.split('=', 1) for line in (run / 'resolved_run.cfg').read_text().splitlines())
        expected = DRIVER.solver_values(resolved, run, run / 'initial_particles.csv', 2)
        self.assertAlmostEqual(float(values['adhesion_work']), expected['adhesion_work'])
        self.assertLess(float(values['adhesion_work']), config['interaction']['adhesion_work_J_m2'])
        self.assertEqual(values['max_steps'], '2')
        manifest = json.loads((run / 'manifest.json').read_text())
        self.assertEqual(manifest['config'], resolved)
        self.assertEqual(manifest['ranks'], 1)
        self.assertTrue((run / 'summarize_particle_solver.py').is_file())
        self.assertTrue((run / 'particles_to_paraview.py').is_file())

    @unittest.skipUnless(shutil.which('c++'), 'C++ compiler unavailable')
    def test_actual_restart_uses_saved_cmc_in_place_without_replacing_snapshot(self):
        config=json.loads((ROOT/'slurry/cases/gr_CMC.json').read_text())
        config['cmc']['net_potential']['contact_force_N']*=.37
        output=self.root/'saved_run';run=output/'gr_cmc/g007_10'
        saved=output/'input/cases/gr_cmc.json';saved.parent.mkdir(parents=True)
        COMMON.write_json(saved,config)
        stale_driver=output/'input/drivers/gr_re2/run_graphite.py'
        stale_driver.parent.mkdir(parents=True)
        stale_driver.write_text('raise RuntimeError("Must not execute the old snapshotted driver")\n')
        point=checkpoint(run,config,8,10)
        data=json.loads((point/'checkpoint.json').read_text())
        data['output_bytes']={'history.csv':9,'particles.csv':9}
        COMMON.write_json(point/'checkpoint.json',data)
        for name in ('history.csv','particles.csv'):(run/name).write_text('header\n8\n12\n')
        (run/'initial_particles.csv').write_bytes(b'x')
        (run/'solver.log').write_text('old solver output\n')
        COMMON.write_json(run/'effective_config.json',DRIVER.resolve(config,10)[0])
        COMMON.write_json(run/'manifest.json',{'saved_run':True})
        before={str(path.relative_to(output/'input')):path.read_bytes()
                for path in (output/'input').rglob('*') if path.is_file()}
        executable=self.root/'fake_restart_slurry'
        info={'mpi_enabled':False,'rough_contact':True,'local_gap_adhesion':True,
              'surface_adhesion_version':1,'pass_max_version':1,'cmc_contact_version':1,
              'free_cmc_repulsion_version':1,'free_cmc_inner_repulsion_version':1,
              'cmc_net_potential_version':1,'cmc_coordination_version':1,
              'current_adhesion_rolling_version':1,'pure_gr_checkpoint_version':1,
              'particle_solver':'petsc','engines':['pure_gr','gr_cmc']}
        executable.write_text(
            '#!'+sys.executable+'\n'
            'import json, os, pathlib, shutil, sys\n'
            'if sys.argv[1:] == ["--build-info"]:\n'
            '    print('+repr(json.dumps(info))+')\n'
            'else:\n'
            '    root=pathlib.Path.cwd()\n'
            '    cfg=dict(line.split("=",1) for line in pathlib.Path(sys.argv[2]).read_text().splitlines())\n'
            '    observed={"cwd":str(root),"engine":os.environ.get("SLURRY_ENGINE"),"config":cfg,\n'
            '              "history_before":(root/"history.csv").read_text()}\n'
            '    (root/"observed_restart.json").write_text(json.dumps(observed))\n'
            '    old=pathlib.Path(cfg["restart_dir"])\n'
            '    for name in ("history.csv","particles.csv"):\n'
            '        with (root/name).open("a") as stream:stream.write("16\\n")\n'
            '    data=json.loads((old/"checkpoint.json").read_text());data["step"]=16\n'
            '    data["output_bytes"]={name:(root/name).stat().st_size for name in data["output_bytes"]}\n'
            '    new=root/"checkpoints"/"checkpoint_00000000000000000016";new.mkdir()\n'
            '    for name in data["files"]:shutil.copy2(str(old/name),str(new/name))\n'
            '    (new/"checkpoint.json").write_text(json.dumps(data))\n'
            '    (root/"status.json").write_text(json.dumps({"status":"MAX_STEPS"}))\n'
            '    print("resumed solver output")\n')
        executable.chmod(0o755)
        COMMON.write_json(self.root/'build_manifest.json',{
            'build_id':'restart-fixture','build_info':info,
            'executable_sha256':COMMON.digest(executable),'source_inputs':COMMON.source_inputs()})
        result=subprocess.run([sys.executable,str(CONTROLLER),'--restart',str(output),
                               '--executable',str(executable),'--max-steps','16','--end-strain','20'],
                              cwd=str(ROOT),env=self.environment,stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,universal_newlines=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        observed=json.loads((run/'observed_restart.json').read_text())
        self.assertEqual(observed['cwd'],str(run))
        self.assertEqual(observed['engine'],'gr_cmc')
        self.assertEqual(float(observed['config']['cmc_net_contact_force']),
                         config['cmc']['net_potential']['contact_force_N'])
        self.assertEqual(float(observed['config']['end_strain']),20.)
        self.assertEqual(observed['history_before'],'header\n8\n')
        for name in ('history.csv','particles.csv'):
            self.assertEqual((run/name).read_text(),'header\n8\n16\n')
        self.assertIn('old solver output',(run/'solver.log').read_text())
        self.assertIn('resumed solver output',(run/'solver.log').read_text())
        self.assertTrue(point.is_dir())
        self.assertTrue((run/'checkpoints/checkpoint_00000000000000000016/checkpoint.json').is_file())
        after={str(path.relative_to(output/'input')):path.read_bytes()
               for path in (output/'input').rglob('*') if path.is_file()}
        self.assertEqual(before,after,'original input snapshot must remain byte-for-byte intact')
        self.assertFalse((output/'pure_gr').exists())

    def test_cpp_dispatcher_routes_alias_with_generated_registry(self):
        engines = json.loads((ROOT / 'slurry/engines.json').read_text())
        declarations = ['namespace slurry { namespace ' + entry['namespace'] +
                        ' { int runCase(int,char**); } }' for entry in engines]
        entries = ['{"' + entry['id'] + '", &slurry::' + entry['namespace'] + '::runCase}'
                   for entry in engines]
        (self.root / 'slurry_registry.h').write_text(
            '\n'.join(declarations) + '\n'
            'struct SlurryEngine { const char* name; int (*entry)(int,char**); };\n'
            'const SlurryEngine slurry_registry[] = {' + ','.join(entries) + '};\n')
        stubs = self.root / 'stubs.cpp'
        stubs.write_text('#include <iostream>\n' + '\n'.join(
            'namespace slurry { namespace ' + namespace + ' {\n'
            'int runCase(int argc,char** argv) { std::cout << "' + namespace + '";\n'
            'for (int i=1; i<argc; ++i) std::cout << " " << argv[i]; return 0; } } }'
            for namespace in sorted({entry['namespace'] for entry in engines})))
        executable = self.root / 'dispatcher'
        compile_result = subprocess.run(
            [shutil.which('c++'), '-std=c++11', '-I', str(self.root),
             str(ROOT / 'olb-1.9r0/examples/slurry/slurry.cpp'), str(stubs), '-o', str(executable)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        self.assertEqual(compile_result.returncode, 0, compile_result.stderr)
        info = json.loads(subprocess.check_output([str(executable), '--build-info'],
                                                 universal_newlines=True))
        self.assertIn('gr_cmc', info['engines'])
        self.assertEqual(info['pure_gr_checkpoint_version'], 1)
        self.assertEqual(info['free_cmc_repulsion_version'], 1)
        self.assertEqual(info['free_cmc_inner_repulsion_version'], 1)
        dispatched = subprocess.check_output(
            [str(executable), '--engine', 'gr_cmc', '--config', 'fixture.cfg'],
            universal_newlines=True)
        self.assertEqual(dispatched, 'gr_re2 --config fixture.cfg')
        environment_dispatch = subprocess.check_output(
            [str(executable), '--config', 'fixture.cfg'],
            env=dict(self.environment, SLURRY_ENGINE='gr_cmc'), universal_newlines=True)
        self.assertEqual(environment_dispatch, dispatched)


if __name__ == '__main__':
    unittest.main()
