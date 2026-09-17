#!/usr/bin/env python3
"""One build, one binary, all registered cases. Python 3.6+."""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from common import BASE, OLB, digest, fingerprint, git_state, source_inputs, tool_output, write_json

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=['mpi', 'serial'], default='mpi')
    p.add_argument('--jobs', type=int, default=int(os.environ.get('SLURRY_BUILD_JOBS', os.environ.get('SLURM_CPUS_PER_TASK', '4'))))
    a = p.parse_args()
    if a.jobs < 1: p.error('--jobs must be positive')
    if any(c.isspace() for c in str(BASE)):
        p.error('OpenLB Makefiles require an installation path without whitespace')
    compiler = os.environ.get('SLURRY_CXX', 'mpicxx' if a.mode == 'mpi' else 'g++')
    compiler_path = shutil.which(compiler)
    if not compiler_path: p.error('Compiler not found: ' + compiler)
    if not (OLB / 'src/case/case.h').is_file(): p.error('OpenLB 1.9r0 source not found at ' + str(OLB))
    toolchain = {'compiler': compiler_path, 'version': tool_output([compiler_path, '--version']),
                 'mode': a.mode, 'platform': 'CPU_SISD', 'config_sha256': digest(OLB / 'config.mk'),
                 'rules_sha256': digest(OLB / 'rules.mk'), 'makefile_sha256': digest(BASE / 'slurry/Makefile')}
    if a.mode == 'mpi':
        toolchain['mpi'] = tool_output([compiler_path, '--showme'])
    build_root = BASE / 'build/slurry'
    build_root.mkdir(parents=True, exist_ok=True)
    with (build_root / 'build.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        inputs = source_inputs()
        build_id = fingerprint({'source': inputs, 'toolchain': toolchain})
        release = build_root / 'releases' / build_id
        complete = release / 'build_manifest.json'
        if complete.is_file() and (release / 'slurry').is_file():
            old = json.loads(complete.read_text())
            if digest(release / 'slurry') == old['executable_sha256']:
                publish(build_root, release)
                print('BUILD REUSED: ' + str(release / 'slurry')); return 0
        work = build_root / 'work' / fingerprint(toolchain)[:16]
        work.mkdir(parents=True, exist_ok=True)
        engines = json.loads((BASE / 'slurry/engines.json').read_text())
        declarations, entries = [], []
        for e in engines:
            if not re.fullmatch(r'[a-z][a-z0-9_]*', e['id']) or not re.fullmatch(r'[a-z][a-z0-9_]*', e['namespace']):
                raise ValueError('Invalid engine registration')
            declarations.append('namespace slurry { namespace ' + e['namespace'] + ' { int runCase(int,char**); } }')
            entries.append('{"' + e['id'] + '", &slurry::' + e['namespace'] + '::runCase}')
        generated = '\n'.join(declarations) + '\nstruct SlurryEngine { const char* name; int (*entry)(int,char**); };\n'
        generated += 'const SlurryEngine slurry_registry[] = {' + ','.join(entries) + '};\n'
        registry = work / 'slurry_registry.h'
        if not registry.exists() or registry.read_text() != generated: registry.write_text(generated)
        # OpenLB 1.9 has non-inline definitions in aggregate headers. Compile
        # one application translation unit, with case namespaces kept distinct.
        sources = sorted((OLB / 'src/slurry').rglob('*.cpp'))
        sources.append(OLB / 'examples/slurry/slurry.cpp')
        unity_text = '#include <olb.h>\n' + ''.join('#include ' + json.dumps(str(s)) + '\n' for s in sources)
        unity = work / 'slurry_application.cpp'
        if not unity.exists() or unity.read_text() != unity_text: unity.write_text(unity_text)
        # Native embedded dependencies; no upstream source or config is edited.
        subprocess.run(['make', '-C', str(OLB / 'external'), '-j' + str(a.jobs), 'lib', 'zlib', 'tinyxml2',
                        'CXX=' + compiler_path, 'CC=gcc'], check=True)
        # A cache miss rebuilds compiled sources even if an archive preserved an
        # older mtime. Exact-content cache hits above need no recompilation.
        command = ['make', '--always-make', '-f', str(BASE / 'slurry/Makefile'), '-j' + str(a.jobs),
                   'OLB_ROOT=' + str(OLB), 'BUILD_DIR=' + str(work), 'CXX=' + compiler_path,
                   'CC=gcc', 'PARALLEL_MODE=' + ('MPI' if a.mode == 'mpi' else 'OFF'),
                   'PLATFORMS=CPU_SISD', 'CUDA_CXX=', 'FEATURES=', 'FLOATING_POINT_TYPE=double']
        subprocess.run(command, cwd=str(BASE), check=True)
        if inputs != source_inputs(): raise RuntimeError('Sources changed while building; retry after edits finish')
        executable = work / 'slurry'
        info = json.loads(tool_output([str(executable), '--build-info']))
        if info['mpi_enabled'] != (a.mode == 'mpi'): raise RuntimeError('MPI build-mode verification failed')
        release.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(executable), str(release / 'slurry'))
        record = {'created_utc': datetime.now(timezone.utc).isoformat(), 'build_id': build_id,
                  'source_inputs': inputs, 'toolchain': toolchain, 'command': command, 'build_info': info,
                  'executable_sha256': digest(release / 'slurry'), 'git': git_state()}
        write_json(complete, record)
        publish(build_root, release)
        print('BUILD COMPLETE: ' + str(build_root / 'current/slurry'))
    return 0

def publish(root, release):
    temporary = root / ('current.tmp.' + str(os.getpid()))
    temporary.symlink_to(release.relative_to(root), target_is_directory=True)
    temporary.replace(root / 'current')

if __name__ == '__main__':
    try: sys.exit(main())
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as e:
        print('BUILD FAILED: ' + str(e), file=sys.stderr); sys.exit(1)
