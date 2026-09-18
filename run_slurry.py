#!/usr/bin/env python3
"""Common batch controller. Original case drivers retain all physical mapping."""
import argparse
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

def positives(text):
    values = [float(x.strip()) for x in text.split(',')]
    if not values or any(not math.isfinite(x) or x <= 0 for x in values):
        raise ValueError('Shear rates must be finite positive numbers')
    return sorted(set(values))

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--settings', type=Path, default=BASE / 'slurry/cases/run.json')
    p.add_argument('--cases', help='pure_cmc,pure_gr,gr_baseline')
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
    p.add_argument('--restart', type=Path, help='Existing gr_baseline checkpoint; same integrated binary/ranks')
    p.add_argument('--output', type=Path)
    p.add_argument('--executable', type=Path, default=BASE / 'build/slurry/current/slurry')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--smoke', action='store_true', help='One CMC rate at 100/s and 20 graphite steps; material/geometry unchanged')
    p.add_argument('--keep-going', action='store_true')
    a = p.parse_args()
    settings = json.loads(a.settings.expanduser().resolve().read_text(encoding='utf-8'))
    cases = [x.strip() for x in a.cases.split(',')] if a.cases else settings['cases']
    registry = {e['id']: e for e in json.loads((BASE / 'slurry/engines.json').read_text())}
    if not cases or len(cases) != len(set(cases)): p.error('Choose at least one case, without duplicates')
    for engine in cases:
        if engine not in registry: p.error('Unknown/unimplemented case: ' + engine + ' (Gr+CMC is not implemented)')
    graphites = [e for e in cases if e != 'pure_cmc']
    if a.config and len(graphites) != 1: p.error('--config requires exactly one graphite case')
    if a.restart and (cases != ['gr_baseline'] or a.shear_rates):
        p.error('--restart requires only gr_baseline and its original configuration/shear rate')
    if any(x is not None for x in (a.target_mach, a.time_step, a.end_strain)) and 'pure_gr' not in cases:
        p.error('--target-mach/--time-step/--end-strain require pure_gr')
    if a.max_steps is not None and a.max_steps < 0: p.error('--max-steps must be nonnegative')
    rates = positives(a.shear_rates) if a.shear_rates else settings.get('shear_rates_s_inv')
    if rates is not None: rates = positives(','.join(str(x) for x in rates))
    if a.smoke and rates is None: rates = [100.0]
    allocation = int(os.environ.get('SLURM_NTASKS', '1'))
    ranks = {}
    for engine in cases:
        value = settings.get('ranks', {}).get(engine, 'allocation')
        value = allocation if value == 'allocation' else int(value)
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
    output = (a.output or default_output).expanduser().resolve()
    if output.exists() and not a.dry_run: p.error('Output already exists: ' + str(output))
    executable = a.executable.expanduser().resolve()
    record = None
    if not a.dry_run:
        record = build_record(executable)
        if max(ranks.values()) > 1 and not record['build_info']['mpi_enabled']: p.error('Multiple ranks require the MPI build')
    configs = {}
    for engine in graphites:
        path = a.config if a.config else a.settings.resolve().parent / settings['graphite_configs'][engine]
        path = path.expanduser().resolve()
        json.loads(path.read_text(encoding='utf-8'))
        configs[engine] = path
    if not a.dry_run:
        output.mkdir(parents=True)
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
            for index, rate in enumerate(rates or [None]):
                name = engine + ('/g%03d_%s' % (index, format(rate, '.12g')) if rate is not None else '')
                result_dir = output / name
                if not a.dry_run:
                    result_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(str(BASE / 'slurry/tools/particles_to_paraview.py'),
                                 str(result_dir / 'particles_to_paraview.py'))
                    if engine == 'pure_gr':
                        shutil.copy2(str(BASE / 'slurry/tools/summarize_particle_solver.py'),
                                     str(result_dir / 'summarize_particle_solver.py'))
                args = prefix + ['--config', str(configs[engine]), '--output', str(result_dir)]
                if rate is not None: args += ['--shear-rate', str(rate)]
                limit = a.max_steps if a.max_steps is not None else (20 if a.smoke else None)
                if limit is not None: args += ['--max-steps' if engine == 'pure_gr' else '--benchmark-steps', str(limit)]
                if engine == 'pure_gr':
                    for key in ('target_mach', 'time_step', 'end_strain'):
                        if getattr(a, key) is not None: args += ['--' + key.replace('_', '-'), str(getattr(a, key))]
                if a.restart: args += ['--restart', str(a.restart.expanduser().resolve())]
                jobs.append({'engine': engine, 'argv': args})
    plan = {'created_utc': datetime.now(timezone.utc).isoformat(), 'jobs': jobs, 'ranks': ranks,
            'output': str(output), 'executable': str(executable), 'git': git_state(),
            'controller_sha256': digest(__file__), 'build_id': record['build_id'] if record else None}
    if a.dry_run:
        print(json.dumps(plan, indent=2)); return 0
    write_json(output / 'batch_manifest.json', plan)
    print('Batch output: ' + str(output), flush=True)
    state = {'status': 'RUNNING', 'results': []}
    write_json(output / 'batch_status.json', state)
    child, current, stop = None, None, []
    def forward(number, frame):
        if number == signal.SIGUSR1 and current != 'gr_baseline':
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
    return 128 + stop[0] if stop else (1 if failed else 0)

if __name__ == '__main__':
    try: sys.exit(main())
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as e:
        print('RUN FAILED: ' + str(e), file=sys.stderr); sys.exit(1)
