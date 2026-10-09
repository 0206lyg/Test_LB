#!/usr/bin/env python3
"""Common batch controller. Original case drivers retain all physical mapping."""
import argparse
import fcntl
import importlib.util
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
from common import BASE, build_record, digest, git_state, write_json

RE2_ENGINES = ('pure_gr', 'gr_cmc')

def positives(text):
    values = [float(x.strip()) for x in text.split(',')]
    if not values or any(not math.isfinite(x) or x <= 0 for x in values):
        raise ValueError('Shear rates must be finite positive numbers')
    return sorted(set(values))

def restart_folder(value):
    path=value.expanduser()
    candidates=[path] if path.is_absolute() else [BASE/path,BASE/'runs'/path]
    for candidate in candidates:
        if candidate.is_dir():
            if candidate.resolve()==(BASE/'runs').resolve():
                raise ValueError('Select one run folder under runs, not the entire runs directory')
            return candidate.resolve()
    raise ValueError('Restart folder not found: '+str(value))

def graphite_restart_driver():
    spec=importlib.util.spec_from_file_location('slurry_gr_restart',BASE/'slurry/drivers/gr_re2/run_graphite.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def saved_json(path):
    return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {}

def restart_context(folder):
    """Locate batch provenance even when the selected path is one checkpoint."""
    for root in (folder,) + tuple(folder.parents):
        if (root/'input/cases').is_dir() or (root/'batch_manifest.json').is_file():
            return root
    # Older/direct runs may have no batch provenance. Keep their directory
    # layout, but require an explicit or saved effective configuration.
    for root in (folder,) + tuple(folder.parents):
        if root.name in RE2_ENGINES:
            return root.parent
    root=folder.parent.parent if (folder/'checkpoint.json').is_file() and folder.parent.name=='checkpoints' else folder
    return root

def acquire_run_lock(root):
    lock=(root/'.run_slurry.lock').open('a')
    try:fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise ValueError('Another batch controller is using '+str(root))
    return lock

def inferred_restart_cases(folder, root, manifest):
    relative=folder.relative_to(root)
    for part in relative.parts:
        if part in RE2_ENGINES:
            return [part]
    found=[]
    for job in manifest.get('jobs', []):
        if job.get('engine') in RE2_ENGINES and job['engine'] not in found:
            found.append(job['engine'])
    for engine in RE2_ENGINES:
        if ((root/engine).is_dir() or any(path.stem.lower()==engine for path in (root/'input/cases').glob('*.json'))):
            if engine not in found:found.append(engine)
    if not found:
        effective=saved_json(root/'effective_config.json')
        if effective:
            found=['gr_cmc' if 'cmc' in effective else 'pure_gr']
    if not found:
        raise ValueError('Cannot infer restart case; retain input/cases and batch_manifest.json or specify --cases and --config')
    return found

def saved_case_config(root, engine, settings):
    candidates=[root/'input/cases'/(engine+'.json')]
    name=settings.get('graphite_configs',{}).get(engine)
    if name:candidates.append(root/'input/cases'/Path(name).name)
    candidates.extend(path for path in (root/'input/cases').glob('*.json') if path.stem.lower()==engine)
    return next((path.resolve() for path in candidates if path.is_file()),None)

def argument_value(argv, option):
    if option not in argv:return None
    index=argv.index(option)+1
    if index==len(argv):raise ValueError('Saved job has no value for '+option)
    return argv[index]

def saved_restart_job(manifest, engine, directory, root):
    archives=sorted((root/'restart_attempts').glob('*/batch_manifest.json'),reverse=True)
    for source in [manifest]+archives:
        previous=saved_json(source) if isinstance(source,Path) else source
        for job in previous.get('jobs', []):
            if job.get('engine')!=engine:continue
            output=argument_value(job.get('argv',[]),'--output')
            if not output:continue
            old=Path(output)
            if old.resolve()==directory:return job
            try:
                if old.relative_to(Path(previous['output']))==directory.relative_to(root):return job
            except (KeyError,ValueError):pass
    return {}

def restart_numerics(args, current, job):
    """Retain original numerical CLI overrides, never an old strain endpoint."""
    values={key:getattr(args,key) for key in ('target_mach','time_step')}
    previous=job.get('input_numerics')
    for key,config_key in (('target_mach','target_mach'),('time_step','time_step_s')):
        original=argument_value(job.get('argv',[]),'--'+key.replace('_','-'))
        if values[key] is None and original is not None:
            if previous is None or current.get('numerics',{}).get(config_key)==previous.get(config_key):
                values[key]=float(original)
    return values

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--settings', type=Path, default=BASE / 'slurry/cases/run.json')
    p.add_argument('--cases', help='pure_cmc,pure_gr,gr_cmc,gr_baseline')
    p.add_argument('--shear-rates', help='Comma-separated s^-1; omitted means each existing case keeps its defaults')
    p.add_argument('--concentrations')
    p.add_argument('--resolution', type=int)
    p.add_argument('--points', type=int)
    p.add_argument('--ranks', type=int, help='Override ranks for every selected case')
    p.add_argument('--cmc-ranks', type=int)
    p.add_argument('--config', type=Path, help='Graphite JSON; requires exactly one selected graphite case')
    p.add_argument('--max-steps', type=int, help='Graphite step cap; CMC convergence settings are retained')
    p.add_argument('--cmc-max-steps', type=int)
    p.add_argument('--target-mach', type=float, help='Existing RE2 option')
    p.add_argument('--time-step', type=float, help='Existing RE2 option in seconds')
    p.add_argument('--end-strain', type=float, help='Existing RE2 option')
    p.add_argument('--restart', type=Path, help='Previous run name or path; infer saved graphite cases and continue in their original directories')
    p.add_argument('--output', type=Path)
    p.add_argument('--executable', type=Path, default=BASE / 'build/slurry/current/slurry')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--smoke', action='store_true', help='One CMC rate at 100/s and 20 graphite steps; material/geometry unchanged')
    p.add_argument('--keep-going', action='store_true')
    a = p.parse_args()
    settings={};restart_root=None;saved_settings={};saved_manifest={};run_lock=None
    if a.restart:
        a.restart=restart_folder(a.restart)
        restart_root=restart_context(a.restart)
        if not a.dry_run and a.cases!='gr_baseline':run_lock=acquire_run_lock(restart_root)
        saved_settings=saved_json(restart_root/'input/run_settings.json')
        saved_manifest=saved_json(restart_root/'batch_manifest.json')
        if saved_settings:settings=saved_settings
    if not a.restart or (a.cases=='gr_baseline' and not saved_settings):
        settings=json.loads(a.settings.expanduser().resolve().read_text(encoding='utf-8'))
    cases = [x.strip() for x in a.cases.split(',')] if a.cases else (inferred_restart_cases(a.restart,restart_root,saved_manifest) if a.restart else settings['cases'])
    registry = {e['id']: e for e in json.loads((BASE / 'slurry/engines.json').read_text())}
    if not cases or len(cases) != len(set(cases)): p.error('Choose at least one case, without duplicates')
    for engine in cases:
        if engine not in registry: p.error('Unknown case: ' + engine)
    if a.restart:
        selected_engine=next((part for part in a.restart.relative_to(restart_root).parts if part in RE2_ENGINES),None)
        if selected_engine and cases!=[selected_engine]:
            p.error('Selected restart directory belongs to '+selected_engine+'; requested cases do not match')
    graphites = [e for e in cases if e != 'pure_cmc']
    if a.config and len(graphites) != 1: p.error('--config requires exactly one graphite case')
    if a.restart and (not all(engine in RE2_ENGINES for engine in cases) and cases!=['gr_baseline'] or a.shear_rates or a.smoke):
        p.error('--restart requires saved pure_gr/gr_cmc cases or one explicit gr_baseline case; it retains saved shear rates and cannot use --smoke')
    in_place=bool(a.restart and all(engine in RE2_ENGINES for engine in cases))
    restart_runs={};restart_skipped=[];restart_unstarted=[];restart_driver=None
    if in_place:
        restart_driver=graphite_restart_driver()
        for engine in cases:
            # A mixed batch may stop before starting one of its saved cases.
            # Skip only that absence, never a malformed completed checkpoint.
            scope=a.restart/engine if any((a.restart/name).is_dir() for name in RE2_ENGINES) else a.restart
            complete=(scope/'checkpoint.json').is_file() or any(
                path.parent.name.startswith('checkpoint_') and path.parent.name[11:].isdigit()
                and path.parent.parent.name=='checkpoints'
                and not any(part in ('restart_attempts','restart_archive') for part in path.relative_to(scope).parts)
                for path in scope.rglob('checkpoint.json'))
            if not complete and not a.cases:
                restart_unstarted.append(str(scope));continue
            restart_runs[engine]=restart_driver.restart_targets(a.restart,engine)
        cases=[engine for engine in cases if engine in restart_runs]
        graphites=list(cases)
        if not cases:p.error('No completed graphite checkpoints in '+str(a.restart))
    if any(x is not None for x in (a.target_mach, a.time_step, a.end_strain)) and not any(e in RE2_ENGINES for e in cases):
        p.error('--target-mach/--time-step/--end-strain require pure_gr or gr_cmc')
    if a.max_steps is not None and a.max_steps < 0: p.error('--max-steps must be nonnegative')
    rates = positives(a.shear_rates) if a.shear_rates else settings.get('shear_rates_s_inv')
    if rates is not None: rates = positives(','.join(str(x) for x in rates))
    if a.restart:rates=None
    if a.smoke and rates is None: rates = [100.0]
    allocation = int(os.environ.get('SLURM_NTASKS', '1'))
    ranks = {}
    for engine in cases:
        value = settings.get('ranks', {}).get(engine, 'allocation')
        value = allocation if value == 'allocation' else int(value)
        if in_place:
            saved_ranks={checkpoint['ranks'] for _,checkpoint in restart_runs[engine]}
            if len(saved_ranks)!=1:p.error('Saved '+engine+' checkpoints use different MPI ranks; restart each case directory separately')
            value=saved_ranks.pop()
        if a.ranks is not None: value = a.ranks
        if engine == 'pure_cmc' and a.cmc_ranks is not None: value = a.cmc_ranks
        if value < 1 or ('SLURM_NTASKS' in os.environ and value > allocation): p.error('Invalid MPI ranks for ' + engine)
        ranks[engine] = value
    cmc = dict(settings.get('cmc', {}))
    for name in ('concentrations', 'resolution', 'points'):
        if getattr(a, name) is not None: cmc[name] = getattr(a, name)
    if 'pure_cmc' in cases and ranks['pure_cmc'] > int(cmc.get('resolution', 64)) // 4:
        p.error('CMC requires ranks <= resolution/4; use --cmc-ranks 1')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    default_output = BASE / 'runs' / ('slurry_' + os.environ.get('SLURM_JOB_ID', 'local') + '_' + stamp)
    output = (restart_root if in_place else (a.output or default_output)).expanduser().resolve()
    if in_place and a.output and a.output.expanduser().resolve()!=output:
        p.error('--restart continues in '+str(output)+'; omit --output or use that same directory')
    if output.exists() and not a.dry_run and not in_place: p.error('Output already exists: ' + str(output))
    executable = a.executable.expanduser().resolve()
    record = None
    if not a.dry_run:
        record = build_record(executable)
        if max(ranks.values()) > 1 and not record['build_info']['mpi_enabled']: p.error('Multiple ranks require the MPI build')
    configs = {}
    for engine in graphites:
        path = a.config if a.config else (saved_case_config(restart_root,engine,saved_settings) if in_place
                                         else a.settings.resolve().parent / settings['graphite_configs'][engine])
        if path:
            path = path.expanduser().resolve()
            json.loads(path.read_text(encoding='utf-8'))
        configs[engine] = path
    if in_place:
        for engine,targets in restart_runs.items():
            pending=[]
            for checkpoint_dir,checkpoint in targets:
                if checkpoint_dir.parent.name!='checkpoints':
                    raise ValueError('In-place batch restart requires the original run/checkpoints directory: '+str(checkpoint_dir))
                directory=checkpoint_dir.parent.parent
                if restart_root!=directory and restart_root not in directory.parents:
                    raise ValueError('Checkpoint is outside the selected run: '+str(checkpoint_dir))
                path=configs[engine] or directory/'effective_config.json'
                if not path.is_file():
                    raise ValueError('No saved configuration for '+engine+' in '+str(restart_root)+'; retain input/cases or supply --config')
                current=json.loads(path.read_text(encoding='utf-8'))
                previous={} if a.config else saved_restart_job(saved_manifest,engine,directory,restart_root)
                options=restart_numerics(a,current,previous)
                cfg,meta=restart_driver.resolve(current,checkpoint['shear_rate_s_inv'],a.max_steps or 0,
                                               options['target_mach'],options['time_step'],a.end_strain)
                if restart_driver.validate_restart(cfg,meta,checkpoint,ranks[engine],a.max_steps or 0,allow_complete=True):
                    pending.append({'checkpoint':checkpoint_dir,'directory':directory,'rate':checkpoint['shear_rate_s_inv'],
                                    'config':path,'numerics':options,'input_numerics':current.get('numerics',{})})
                else:restart_skipped.append(str(checkpoint_dir))
            restart_runs[engine]=pending
        if not any(restart_runs.values()):
            p.error('All selected checkpoints have reached the requested endpoint; increase flow.end_strain or --max-steps')
    if not a.dry_run and not in_place:
        output.mkdir(parents=True)
        run_lock=acquire_run_lock(output)
        snapshot = output / 'input'
        shutil.copytree(str(BASE / 'slurry/drivers'), str(snapshot / 'drivers'), ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        (snapshot / 'cases').mkdir()
        table = BASE / 'slurry/cases/cmc250k_cross_parameters.csv'
        shutil.copy2(str(table), str(snapshot / 'cases' / table.name))
        for engine, path in list(configs.items()):
            destination = snapshot / 'cases' / (engine + '.json')
            shutil.copy2(str(path), str(destination)); configs[engine] = destination
        drivers = snapshot / 'drivers'
        parameters = snapshot / 'cases/cmc250k_cross_parameters.csv'
        write_json(snapshot / 'build_manifest.json', record)
        write_json(snapshot / 'run_settings.json', settings)
    else:
        drivers = BASE / 'slurry/drivers'
        parameters = BASE / 'slurry/cases/cmc250k_cross_parameters.csv'
    jobs = []
    for engine in cases:
        entry = registry[engine]
        driver = drivers / Path(entry['driver']).relative_to('slurry/drivers')
        prefix = [sys.executable, str(driver), '--executable', str(executable), '--ranks', str(ranks[engine])]
        if engine == 'pure_cmc':
            args = prefix + ['--parameters', str(parameters), '--output', str(output / engine)]
            for key, value in cmc.items(): args += ['--' + key.replace('_', '-'), str(value)]
            if rates: args += ['--shear-rates', ','.join(str(x) for x in rates)]
            if a.cmc_max_steps: args += ['--max-steps', str(a.cmc_max_steps)]
            jobs.append({'engine': engine, 'argv': args})
        else:
            selections=restart_runs[engine] if in_place else [{'rate':rate} for rate in (rates or [None])]
            for index, selected in enumerate(selections):
                rate=selected['rate'];restart_case=selected.get('checkpoint')
                name = engine + ('/g%03d_%s' % (index, format(rate, '.12g')) if rate is not None else '')
                result_dir = selected['directory'] if in_place else output / name
                if not a.dry_run and not in_place:
                    result_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(str(BASE / 'slurry/tools/particles_to_paraview.py'),
                                 str(result_dir / 'particles_to_paraview.py'))
                    if engine in RE2_ENGINES:
                        shutil.copy2(str(BASE / 'slurry/tools/summarize_particle_solver.py'),
                                     str(result_dir / 'summarize_particle_solver.py'))
                config=selected['config'] if in_place else configs[engine]
                args = prefix + ['--config', str(config), '--output', str(result_dir)]
                if rate is not None: args += ['--shear-rate', str(rate)]
                limit = a.max_steps if a.max_steps is not None else (20 if a.smoke else None)
                if limit is not None: args += ['--max-steps' if engine in RE2_ENGINES else '--benchmark-steps', str(limit)]
                if engine in RE2_ENGINES:
                    for key in ('target_mach', 'time_step', 'end_strain'):
                        value=selected['numerics'][key] if in_place and key!='end_strain' else getattr(a,key)
                        if value is not None: args += ['--' + key.replace('_', '-'), str(value)]
                if a.restart: args += ['--restart', str(restart_case or a.restart)]
                jobs.append({'engine': engine, 'argv': args,
                             'input_numerics':json.loads(config.read_text(encoding='utf-8')).get('numerics',{})})
    plan = {'created_utc': datetime.now(timezone.utc).isoformat(), 'jobs': jobs, 'ranks': ranks,
            'output': str(output), 'executable': str(executable), 'git': git_state(),
            'restart_source':str(a.restart) if a.restart else None,
            'restart_skipped_complete':restart_skipped,
            'restart_skipped_unstarted':restart_unstarted,
            'controller_sha256': digest(__file__), 'build_id': record['build_id'] if record else None}
    if a.dry_run:
        print(json.dumps(plan, indent=2)); return 0
    if in_place and (output/'batch_manifest.json').is_file():
        archive=output/'restart_attempts'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
        archive.mkdir(parents=True)
        for name in ('batch_manifest.json','batch_status.json'):
            if (output/name).is_file():shutil.copy2(str(output/name),str(archive/name))
    write_json(output / 'batch_manifest.json', plan)
    print('Batch output: ' + str(output), flush=True)
    state = {'status': 'RUNNING', 'results': []}
    write_json(output / 'batch_status.json', state)
    child, current, stop = None, None, []
    def forward(number, frame):
        if number == signal.SIGUSR1 and current not in RE2_ENGINES + ('gr_baseline',):
            print('Slurm time notice received; selected case retains its original stop behavior.', flush=True)
            return
        stop.append(number)
        if child is not None and child.poll() is None: child.send_signal(number)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1): signal.signal(sig, forward)
    for job in jobs:
        current = job['engine']
        environment = dict(os.environ, SLURRY_ENGINE=current)
        print('CASE: ' + current, flush=True)
        child = subprocess.Popen(job['argv'], env=environment)
        code = child.wait()
        state['results'].append({'engine': current, 'exit_code': code, 'argv': job['argv']})
        write_json(output / 'batch_status.json', state)
        if stop or (code != 0 and not a.keep_going): break
    failed = any(r['exit_code'] != 0 for r in state['results'])
    state['status'] = 'STOPPED' if stop else ('FAILED' if failed else 'COMPLETED')
    write_json(output / 'batch_status.json', state)
    print('BATCH ' + state['status'] + ': ' + str(output), flush=True)
    return (0 if stop[0]==signal.SIGUSR1 and not failed else 128+stop[0]) if stop else (1 if failed else 0)

if __name__ == '__main__':
    try: sys.exit(main())
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as e:
        print('RUN FAILED: ' + str(e), file=sys.stderr); sys.exit(1)
