"""Same-directory checkpoint rollback and continuation without a bulk simulation."""
import copy
import fcntl
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
DRIVER=ROOT/'slurry/drivers/gr_re2/run_graphite.py'
SPEC=importlib.util.spec_from_file_location('inplace_graphite',DRIVER)
RUNNER=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(RUNNER)
FOOTER='</Collection>\n</VTKFile>\n'
BUILD={'mpi_enabled':False,'rough_contact':True,'particle_solver':'petsc',
       'pass_max_version':1,'surface_adhesion_version':1,'pure_gr_checkpoint_version':1,
       'free_cmc_repulsion_version':1,'free_cmc_inner_repulsion_version':1,
       'cmc_net_potential_version':1,'cmc_coordination_version':1,
       'current_adhesion_rolling_version':1,'vtk_restart_version':1}


def case():
    value=json.loads((ROOT/'slurry/cases/gr_CMC.json').read_text())
    value['particles']['count']=1
    value['output'].update(sample_every_steps=2,vtk_every_steps=4)
    return value


def save_checkpoint(run,step,cfg,prefix):
    resolved,meta=RUNNER.resolve(cfg)
    values=RUNNER.solver_values(resolved,run,run/'initial_particles.csv',0)
    immutable={}
    for physical,version,keys in (
            ('surface_adhesion','interaction_model_version',RUNNER.SURFACE_CHECKPOINT_KEYS),
            ('free_cmc_repulsion_pressure','free_cmc_repulsion_version',RUNNER.FREE_CMC_CHECKPOINT_KEYS),
            ('free_cmc_inner_repulsion_work','free_cmc_inner_repulsion_version',RUNNER.INNER_CMC_CHECKPOINT_KEYS),
            ('cmc_net_blend','cmc_net_potential_version',RUNNER.NET_CMC_CHECKPOINT_KEYS),
            ('cmc_coordination_enabled','cmc_coordination_version',RUNNER.CMC_COORDINATION_CHECKPOINT_KEYS),
            ('current_adhesion_rolling','current_adhesion_rolling_version',RUNNER.CURRENT_ADHESION_ROLLING_CHECKPOINT_KEYS)):
        if values.get(physical,0)>0:
            immutable.update({key:values[key] for key in keys if key!=version})
            immutable[version]=1
    immutable.update(shear_rate=100,dt_s=meta['dt_s'],particle_count=1,ranks=1)
    path=run/'checkpoints'/('checkpoint_%020d'%step);path.mkdir(parents=True)
    contents={'state.bin':b'old-state-'+str(step).encode(),'lattice_rank_0.bin':b'lattice',
              'initial_particles.csv':b'original particles\n'}
    for name,value in contents.items():(path/name).write_bytes(value)
    data={'engine':'pure_gr','format_version':1,'complete':True,'step':step,'ranks':1,
          'dt_s':meta['dt_s'],'shear_rate_s_inv':100,'immutable_config':immutable,
          'files':{name:len(value) for name,value in contents.items()},
          'output_bytes':{name:len(prefix) for name in ('history.csv','particles.csv')}}
    RUNNER.write_json(path/'checkpoint.json',data)
    return path,data


class InPlaceRestartTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='graphite-inplace-')
        self.root=Path(self.temp.name)
        self.run=self.root/'gr_cmc'/'g000_100';self.run.mkdir(parents=True)
        self.cfg=case()
        self.prefix=b'step,value\n0,initial\n4,saved\n8,saved\n'
        self.extra=b'12,uncommitted\n'
        self.checkpoint,self.data=save_checkpoint(self.run,8,self.cfg,self.prefix)
        for name in ('history.csv','particles.csv'):(self.run/name).write_bytes(self.prefix+self.extra)
        (self.run/'initial_particles.csv').write_bytes(b'original particles\n')
        RUNNER.write_json(self.run/'effective_config.json',self.cfg)
        for name in ('manifest.json','mapping.json','status.json','driver_status.json'):
            RUNNER.write_json(self.run/name,{'previous':True})
        (self.run/'resolved_run.cfg').write_text('previous=1\n')
        (self.run/'solver.log').write_text('previous solver output\n')
        (self.run/'STOP_REQUEST').touch()
        (self.run/'keep-user-note.txt').write_text('keep this')
        self.vtk=self.run/'vtk'/'vtkData';(self.vtk/'data').mkdir(parents=True)
        entries=[]
        for frame in (1,2,3):
            stem='graphite_iT%07d'%frame
            (self.vtk/'data'/(stem+'.vtm')).write_text('frame '+str(frame))
            (self.vtk/'data'/(stem+'iC00000.vti')).write_text('block '+str(frame))
            entries.append('<DataSet timestep="%d" group="" part="" file="data/%s.vtm"/>'%(frame,stem))
        (self.vtk/'graphite.pvd').write_text(
            '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n<Collection>\n'
            +'\n'.join(entries)+'\n'+FOOTER)
        self.executable=self.root/'fake_graphite'
        self.fake_executable()

    def tearDown(self):
        self.temp.cleanup()

    def fake_executable(self,build=None):
        # Exercise the real driver and its filesystem behavior. This executable
        # only models C++ append/publication semantics, not particle dynamics.
        self.executable.write_text('#!'+sys.executable+'\n'+'''import json,pathlib,shutil,sys
BUILD='''+repr(BUILD if build is None else build)+'''
if sys.argv[1:]==['--build-info']:
 print(json.dumps(BUILD));raise SystemExit(0)
values=dict(line.split('=',1) for line in pathlib.Path(sys.argv[2]).read_text().splitlines())
out=pathlib.Path(values['output_dir']);cp=pathlib.Path(values['restart_dir'])
saved=json.loads((cp/'checkpoint.json').read_text());start=saved['step'];stop=int(values['max_steps'])
assert not (out/'STOP_REQUEST').exists()
assert (out/'initial_particles.csv').read_bytes()==(cp/'initial_particles.csv').read_bytes()
for step in range(start+1,stop+1):
 if step%int(values['sample_every'])==0 or step==stop:
  for name in ('history.csv','particles.csv'):
   with (out/name).open('a') as stream:stream.write(str(step)+',continued\\n')
 period=int(values['vtk_every'])
 if period and step%period==0:
  folder=out/'vtk'/'vtkData';master=folder/'graphite.pvd'
  text=master.read_text();assert 'vtk_time_index="simulation_step"' in text
  stem='graphite_iT%07d'%step
  (folder/'data'/(stem+'.vtm')).write_text('new frame')
  (folder/'data'/(stem+'iC00000.vti')).write_text('new block')
  footer='</Collection>\\n</VTKFile>\\n'
  entry='<DataSet timestep="%d" group="" part="" file="data/%s.vtm"/>\\n'%(step,stem)
  assert text.endswith(footer);master.write_text(text[:-len(footer)]+entry+footer)
new=out/'checkpoints'/('checkpoint_%020d'%stop);new.mkdir()
for name in saved['files']:shutil.copy2(cp/name,new/name)
saved['step']=stop;saved['output_bytes']={name:(out/name).stat().st_size for name in ('history.csv','particles.csv')}
(new/'checkpoint.json').write_text(json.dumps(saved))
(out/'latest_checkpoint.txt').write_text('checkpoints/'+new.name+'\\n')
(out/'status.json').write_text(json.dumps({'status':'MAX_STEPS','step':stop}))
print('continued through step '+str(stop))
''')
        self.executable.chmod(0o755)

    def launch(self,*args):
        return subprocess.run([sys.executable,str(DRIVER),'--restart',str(self.run),
                               '--executable',str(self.executable),'--max-steps','12']+list(args),
                              stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)

    def test_same_file_history_keeps_exact_committed_prefix(self):
        self.assertTrue(RUNNER.copy_history_prefix(self.checkpoint,self.data,self.run))
        for name in ('history.csv','particles.csv'):
            self.assertEqual((self.run/name).read_bytes(),self.prefix)
        self.assertEqual((self.checkpoint/'state.bin').read_bytes(),b'old-state-8')

    def test_rollback_archives_tails_future_checkpoints_and_stop_request(self):
        future,_=save_checkpoint(self.run,12,self.cfg,self.prefix+self.extra)
        partial=self.run/'checkpoints/checkpoint_00000000000000000016.partial';partial.mkdir()
        (self.run/'latest_checkpoint.txt').write_text('checkpoints/'+future.name+'\n')
        plan=RUNNER.inplace_restart_plan(self.checkpoint,self.data,self.run)
        archive=RUNNER.prepare_inplace_restart(self.checkpoint,self.data,self.run,plan)
        for name in ('history.csv','particles.csv'):
            self.assertEqual((self.run/name).read_bytes(),self.prefix)
            self.assertEqual((archive/(name+'.after_checkpoint')).read_bytes(),self.extra)
        self.assertTrue((archive/'STOP_REQUEST').exists())
        self.assertTrue((archive/'status.json').exists())
        self.assertFalse((self.run/'STOP_REQUEST').exists())
        self.assertFalse((self.run/'status.json').exists())
        self.assertTrue((archive/'checkpoints'/future.name).is_dir())
        self.assertTrue((archive/'checkpoints'/partial.name).is_dir())
        self.assertFalse(future.exists());self.assertFalse(partial.exists())
        self.assertEqual((self.run/'latest_checkpoint.txt').read_text(),'checkpoints/'+self.checkpoint.name+'\n')
        self.assertEqual(RUNNER.latest_checkpoint(self.run)[0],self.checkpoint)
        self.assertEqual([path for path,_ in RUNNER.restart_targets(self.root,'gr_cmc')],[self.checkpoint])
        self.assertEqual((self.run/'solver.log').read_text(),'previous solver output\n')
        self.assertEqual((self.run/'keep-user-note.txt').read_text(),'keep this')
        self.assertEqual((self.run/'initial_particles.csv').read_bytes(),b'original particles\n')
        root=ET.parse(self.vtk/'graphite.pvd').getroot()
        self.assertEqual(root.get('vtk_time_index'),'simulation_step')
        self.assertEqual([int(item.get('timestep')) for item in root.find('Collection')],[4,8])
        self.assertTrue((self.vtk/'graphite.pvd').read_text().endswith(FOOTER))
        self.assertTrue((self.vtk/'data/graphite_iT0000002.vtm').exists())
        self.assertFalse((self.vtk/'data/graphite_iT0000003.vtm').exists())
        self.assertTrue((archive/'vtk/vtkData/data/graphite_iT0000003iC00000.vti').exists())

    def test_repeated_direct_restart_reuses_saved_config_and_continuous_outputs(self):
        original=(self.run/'initial_particles.csv').read_bytes()
        # Changing VTK cadence is allowed; old file indices must not be reinterpreted.
        override=copy.deepcopy(self.cfg);override['output']['vtk_every_steps']=3
        path=self.root/'override.json';RUNNER.write_json(path,override)
        result=self.launch('--config',str(path))
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        result=self.launch('--max-steps','16')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        for name in ('history.csv','particles.csv'):
            steps=[int(line.split(',')[0]) for line in (self.run/name).read_text().splitlines()[1:]]
            self.assertEqual(steps,[0,4,8,10,12,14,16])
        master=ET.parse(self.vtk/'graphite.pvd').getroot()
        self.assertEqual([int(item.get('timestep')) for item in master.find('Collection')],[4,8,9,12,15])
        self.assertEqual((self.run/'initial_particles.csv').read_bytes(),original)
        log=(self.run/'solver.log').read_text()
        self.assertIn('previous solver output',log)
        self.assertIn('continued through step 12',log);self.assertIn('continued through step 16',log)
        manifest=json.loads((self.run/'manifest.json').read_text())
        self.assertTrue(manifest['restart_in_place']);self.assertEqual(manifest['restart_step'],12)
        self.assertEqual(manifest['config']['output']['vtk_every_steps'],3)
        self.assertEqual(len(list((self.run/'restart_attempts').iterdir())),2)

    def test_preflight_rejects_changed_physics_before_touching_outputs(self):
        changed=copy.deepcopy(self.cfg);changed['cmc']['net_potential']['contact_force_N']*=2
        path=self.root/'changed.json';RUNNER.write_json(path,changed)
        result=self.launch('--config',str(path))
        self.assertNotEqual(result.returncode,0)
        self.assertIn('physical settings',result.stderr)
        self.assertEqual((self.run/'history.csv').read_bytes(),self.prefix+self.extra)
        self.assertTrue((self.run/'STOP_REQUEST').exists())
        self.assertFalse((self.run/'restart_attempts').exists())

    def test_unsupported_vtk_binary_and_active_rate_lock_leave_run_intact(self):
        old=dict(BUILD);old.pop('vtk_restart_version');self.fake_executable(old)
        result=self.launch()
        self.assertNotEqual(result.returncode,0);self.assertIn('vtk_restart_version',result.stderr)
        self.fake_executable()
        with (self.run/'.run_graphite.lock').open('a+') as lock:
            fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
            result=self.launch()
        self.assertNotEqual(result.returncode,0);self.assertIn('already running',result.stderr)
        self.assertEqual((self.run/'history.csv').read_bytes(),self.prefix+self.extra)
        self.assertFalse((self.run/'restart_attempts').exists())

    def test_missing_or_changed_saved_files_are_rejected_without_mutation(self):
        (self.run/'initial_particles.csv').write_bytes(b'different')
        with self.assertRaisesRegex(ValueError,'Initial particle file differs'):
            RUNNER.inplace_restart_plan(self.checkpoint,self.data,self.run)
        (self.run/'initial_particles.csv').write_bytes(b'original particles\n')
        (self.run/'particles.csv').write_bytes(b'short')
        with self.assertRaisesRegex(ValueError,'shorter'):
            RUNNER.inplace_restart_plan(self.checkpoint,self.data,self.run)
        (self.run/'particles.csv').unlink()
        with self.assertRaisesRegex(ValueError,'incomplete'):
            RUNNER.inplace_restart_plan(self.checkpoint,self.data,self.run)
        self.assertEqual((self.run/'history.csv').read_bytes(),self.prefix+self.extra)
        self.assertTrue((self.run/'STOP_REQUEST').exists())

    def test_missing_vtk_master_can_rebuild_and_malformed_master_is_not_erased(self):
        master=self.vtk/'graphite.pvd';master.unlink()
        plan=RUNNER.vtk_restart_plan(self.run,8,self.cfg)
        root=ET.fromstring(plan['text'])
        self.assertEqual([int(item.get('timestep')) for item in root.find('Collection')],[4,8])
        master.write_text('<broken')
        with self.assertRaisesRegex(ValueError,'Invalid restart VTK'):
            RUNNER.inplace_restart_plan(self.checkpoint,self.data,self.run)
        self.assertEqual(master.read_text(),'<broken')

    def test_missing_absolute_vtk_master_uses_saved_step_convention(self):
        (self.vtk/'graphite.pvd').unlink()
        for path in list((self.vtk/'data').iterdir()):
            index=int(path.name.split('iT')[1][:7])
            path.rename(path.with_name(path.name.replace('iT%07d'%index,'iT%07d'%(4*index))))
        (self.run/'mapping.json').write_text(json.dumps({'vtk_time_index':'simulation_step'}))
        plan=RUNNER.vtk_restart_plan(self.run,8,self.cfg)
        root=ET.fromstring(plan['text'])
        self.assertEqual([int(item.get('timestep')) for item in root.find('Collection')],[4,8])
        self.assertEqual({path.name for path in plan['future']},
                         {'graphite_iT0000012.vtm','graphite_iT0000012iC00000.vti'})

    def test_missing_master_after_legacy_migration_refuses_ambiguous_times(self):
        plan=RUNNER.inplace_restart_plan(self.checkpoint,self.data,self.run)
        RUNNER.prepare_inplace_restart(self.checkpoint,self.data,self.run,plan)
        (self.vtk/'graphite.pvd').unlink()
        (self.run/'mapping.json').write_text(json.dumps({'vtk_time_index':'simulation_step'}))
        before={path.name:path.read_bytes() for path in (self.vtk/'data').iterdir()}
        with self.assertRaisesRegex(ValueError,'cannot safely infer'):
            RUNNER.inplace_restart_plan(self.checkpoint,self.data,self.run)
        self.assertEqual({path.name:path.read_bytes() for path in (self.vtk/'data').iterdir()},before)

    def test_dry_run_uses_saved_settings_without_creating_or_changing_files(self):
        before={str(path.relative_to(self.run)):path.read_bytes() for path in self.run.rglob('*') if path.is_file()}
        result=self.launch('--dry-run')
        self.assertEqual(result.returncode,0,result.stderr)
        value=json.loads(result.stdout)
        self.assertEqual(value['config']['cmc']['net_potential']['contact_force_N'],self.cfg['cmc']['net_potential']['contact_force_N'])
        after={str(path.relative_to(self.run)):path.read_bytes() for path in self.run.rglob('*') if path.is_file()}
        self.assertEqual(after,before)


if __name__=='__main__':unittest.main()
