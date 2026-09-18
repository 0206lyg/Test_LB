#!/usr/bin/env python3
"""Compile a standalone C++ contact test with the same PETSc discovery as OpenLB."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
from common import BASE
from petsc_config import discover, preflight


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--compiler', default=os.environ.get('SLURRY_CXX', 'mpicxx'))
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--serial', action='store_true', help='Allow MPIUNI PETSc for standalone checks')
    args, runtime = parser.parse_known_args()
    compiler = shutil.which(args.compiler)
    if not compiler:
        parser.error('Compiler not found: ' + args.compiler)
    source = args.source.resolve()
    if not source.is_file():
        parser.error('Source file not found: ' + str(source))
    target = (args.output or BASE / 'build/petsc-tests' / source.stem).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    config = preflight(discover(), compiler, not args.serial)
    subprocess.run([compiler, '-std=c++20', '-O2', '-Wall', '-DSLURRY_USE_PETSC',
                    '-I' + str(BASE / 'olb-1.9r0/src'),
                    '-I' + str(BASE / 'olb-1.9r0/src/slurry/gr_re2')] + config['cflags'] +
                   [str(source), '-o', str(target)] + config['libs'], check=True)
    print('COMPILED: ' + str(target), flush=True)
    if args.run:
        if runtime and runtime[0] == '--':
            runtime = runtime[1:]
        subprocess.run(config.get('runtime_launcher', []) + [str(target)] + runtime, check=True)
    elif runtime:
        parser.error('Extra arguments require --run')


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        print('TEST BUILD FAILED: ' + str(error), file=sys.stderr)
        sys.exit(1)
