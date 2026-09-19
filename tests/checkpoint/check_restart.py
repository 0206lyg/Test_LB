#!/usr/bin/env python3
"""Compare an uninterrupted contact run with a saved/restored pure_gr run."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[2]
SPEC=importlib.util.spec_from_file_location('graphite_restart_test',ROOT/'slurry/drivers/gr_re2/run_graphite.py')
driver=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(driver)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,default=ROOT/'build/slurry/current/slurry')
    parser.add_argument('--work',type=Path)
    parser.add_argument('--ranks',type=int,default=1)
    args=parser.parse_args()
    work=args.work or Path(tempfile.mkdtemp(prefix='pure-gr-restart-'))
    work.mkdir(parents=True,exist_ok=True)
    config=json.loads((ROOT/'slurry/cases/pure_gr.json').read_text())
    executable=args.executable.resolve()
    build_info=json.loads(subprocess.check_output([str(executable),'--build-info'],universal_newlines=True))
    config['numerics']['particle_solver']=build_info['particle_solver']
    config['particles']['count']=2
    config['geometry']['dx_m']=1.25e-7
    config['flow']['end_strain']=0.01
    config['output'].update(sample_every_steps=1,checkpoint_every_steps=4,
                            checkpoint_every_seconds=0,checkpoint_keep=2)
    config_path=work/'config.json';driver.write_json(config_path,config)
    cfg,_=driver.resolve(config)
    environment=dict(os.environ,SLURRY_ENGINE='pure_gr')
    prefix=[shutil.which('mpirun'),'-np',str(args.ranks)] if args.ranks>1 else []
    for name,steps in [('whole',16),('part',8)]:
        output=work/name;output.mkdir()
        particles=output/'initial_particles.csv'
        particles.write_text('id,x_m,y_m,z_m,angle_x_deg,angle_y_deg,angle_z_deg\n'
                             '0,5e-6,4.799e-6,5e-6,90,0,0\n'
                             '1,5e-6,5.201e-6,5e-6,90,0,0\n')
        values=driver.solver_values(cfg,output,particles,steps)
        values_path=output/'run.cfg'
        values_path.write_text(''.join('%s=%s\n'%item for item in values.items()))
        with (output/'run.log').open('w') as log:
            subprocess.run(prefix+[str(executable),'--engine','pure_gr','--config',str(values_path)],
                           stdout=log,stderr=subprocess.STDOUT,env=environment,check=True)
    # Exercise the same common controller and driver as run_slurry_cpu.sbatch.
    with (work/'restart.log').open('w') as log:
        subprocess.run([os.sys.executable,str(ROOT/'slurry/tools/run_slurry.py'),
                        '--restart',str(work/'part'),'--config',str(config_path),
                        '--ranks',str(args.ranks),'--max-steps','16','--output',str(work/'resumed'),
                        '--executable',str(executable)],stdout=log,stderr=subprocess.STDOUT,check=True)
    resumed=work/'resumed/pure_gr/g000_100'
    timing={'wall_seconds','steps_per_second','fluid_seconds','map_seconds',
            'coupling_seconds','particle_seconds','output_seconds'}
    for name in ('history.csv','particles.csv'):
        with (work/'whole'/name).open() as f:whole=list(csv.DictReader(f))
        with (resumed/name).open() as f:continued=list(csv.DictReader(f))
        assert len(whole)==len(continued),(name,len(whole),len(continued))
        for i,(a,b) in enumerate(zip(whole,continued)):
            for key in a:
                if key not in timing:
                    assert a[key]==b[key],(name,i,key,a[key],b[key])
        if name=='history.csv':
            assert max(int(row['contact_count']) for row in whole)>0,'Fixture needs active contact history'
            assert max(float(row['contact_elastic_energy_J']) for row in whole)>0,'Fixture needs stored contact elasticity'
    a,_=driver.latest_checkpoint(work/'whole');b,_=driver.latest_checkpoint(resumed)
    for rank in range(args.ranks):
        name='lattice_rank_%d.bin'%rank
        assert (a/name).read_bytes()==(b/name).read_bytes(),'Final fluid distributions differ: '+name
    print('PASS: exact physical CSV rows and final lattice match; active contact history; shared --restart controller; ranks=%d'%args.ranks)
    print('Results:',work)


if __name__=='__main__':main()
