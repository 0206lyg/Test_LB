#!/usr/bin/env python3
"""Exercise SIGUSR1 -> controller -> driver -> collective checkpoint and restore.

Use the part/config.json output of check_restart.py as --source and --config.
Also corrupt one state byte and require the actual solver to reject its checksum.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
SPEC=importlib.util.spec_from_file_location('stop_graphite',ROOT/'slurry/drivers/gr_re2/run_graphite.py')
driver=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(driver)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--executable',type=Path,default=ROOT/'build/slurry/current/slurry')
    parser.add_argument('--work',type=Path)
    args=parser.parse_args()
    work=args.work or Path(tempfile.mkdtemp(prefix='pure-gr-stop-'));work.mkdir(parents=True,exist_ok=True)
    source,data=driver.latest_checkpoint(args.source)
    config=json.loads(args.config.read_text());config['flow']['end_strain']=1
    config_path=work/'config.json';driver.write_json(config_path,config)
    prefix=[sys.executable,str(ROOT/'slurry/tools/run_slurry.py'),'--config',str(config_path),
            '--ranks',str(data['ranks']),'--executable',str(args.executable.resolve())]
    output=work/'stopped'
    result=output/('pure_gr/g000_%s'%format(data['shear_rate_s_inv'],'.12g'))
    with (work/'stop.log').open('w') as log:
        process=subprocess.Popen(prefix+['--restart',str(source),'--max-steps','1000','--output',str(output)],
                                 stdout=log,stderr=subprocess.STDOUT)
        try:
            deadline=time.monotonic()+90
            while not (result/'latest_checkpoint.txt').is_file():
                if process.poll() is not None:raise RuntimeError('Run exited before stop signal; inspect '+str(work/'stop.log'))
                if time.monotonic()>deadline:raise RuntimeError('Timed out waiting for initial restored checkpoint')
                time.sleep(0.05)
            process.send_signal(signal.SIGUSR1)
            assert process.wait(timeout=90)==0,'Controller did not stop cleanly'
        finally:
            if process.poll() is None:process.kill();process.wait()
    assert json.loads((output/'batch_status.json').read_text())['status']=='STOPPED'
    status=json.loads((result/'driver_status.json').read_text())
    assert status['stop_requested'] and status['status']=='CHECKPOINTED',status
    checkpoint,saved=driver.latest_checkpoint(result)
    assert data['step']<saved['step']<1000,saved['step']
    with (work/'continue.log').open('w') as log:
        subprocess.run(prefix+['--restart',str(result),'--max-steps',str(saved['step']+2),
                               '--output',str(work/'continued')],stdout=log,stderr=subprocess.STDOUT,check=True)
    # File size is unchanged, so metadata preflight passes and the C++ checksum
    # must catch corruption before applying any restored state.
    broken=work/'corrupt';broken.mkdir()
    for name in ('checkpoint.json','state.bin','initial_particles.csv'):shutil.copy2(str(checkpoint/name),str(broken/name))
    for rank in range(saved['ranks']):
        name='lattice_rank_%d.bin'%rank;(broken/name).symlink_to(checkpoint/name)
    with (broken/'state.bin').open('r+b') as f:
        f.seek(64);value=f.read(1);f.seek(64);f.write(bytes([value[0]^1]))
    with (work/'corruption.log').open('w') as log:
        failed=subprocess.run(prefix+['--restart',str(broken),'--max-steps',str(saved['step']+2),
                                      '--output',str(work/'rejected')],stdout=log,stderr=subprocess.STDOUT)
    assert failed.returncode!=0 and 'Checkpoint checksum mismatch' in (work/'corruption.log').read_text()
    print('PASS: SIGUSR1 saved at LB boundary, successful subsequent restart, corrupt state rejected')
    print('Results:',work)


if __name__=='__main__':main()
