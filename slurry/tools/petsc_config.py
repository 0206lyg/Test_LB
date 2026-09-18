#!/usr/bin/env python3
"""Discover PETSc and check the exact compiler/MPI combination. Python 3.6+."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile


def output(argv, env=None, cwd=None):
    return subprocess.check_output(argv, env=env, cwd=cwd, stderr=subprocess.STDOUT,
                                   universal_newlines=True).strip()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def wrapper_command(compiler):
    for option in ('--showme', '-show'):
        try:
            value = output([compiler, option])
            if value:
                return value
        except subprocess.CalledProcessError:
            pass
    return output([compiler, '--version'])


def mpi_launcher_prefix(compiler):
    """Use a one-rank launcher from the compiler's MPI installation.

    Direct execution of even a one-rank MPI probe can enter Slurm's PMI
    initialization path. Open MPI's launcher must establish that runtime.
    """
    compiler_directory = Path(compiler).absolute().parent
    for name in ('mpirun', 'mpiexec'):
        candidate = compiler_directory / name
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            return [str(candidate), '-np', '1']
    for name in ('mpirun', 'mpiexec'):
        candidate = shutil.which(name)
        if candidate:
            return [candidate, '-np', '1']
    raise RuntimeError('MPI launcher not found for PETSc runtime check. '
                       'Load the MPI module that supplies the selected compiler: ' + str(compiler))


def discover():
    """Use installed metadata; never guess PETSc's transitive libraries."""
    env = os.environ.copy()
    root = Path(env['PETSC_DIR']).expanduser().resolve() if env.get('PETSC_DIR') else None
    arch = env.get('PETSC_ARCH', '')
    if root:
        prefix = root / arch if arch else root
        paths = [str(prefix / 'lib/pkgconfig'), str(root / 'lib/pkgconfig')]
        env['PKG_CONFIG_PATH'] = os.pathsep.join(paths + [env.get('PKG_CONFIG_PATH', '')])
    pkgconfig = shutil.which('pkg-config')
    package = None
    if pkgconfig:
        for candidate in ('PETSc', 'petsc'):
            try:
                output([pkgconfig, '--exists', candidate], env=env)
                if root:
                    metadata_dir = Path(output([pkgconfig, '--variable=pcfiledir', candidate], env=env)).resolve()
                    if metadata_dir not in [Path(p).resolve() for p in paths]:
                        continue
                package = candidate
                break
            except subprocess.CalledProcessError:
                pass
    if package:
        flags = shlex.split(output([pkgconfig, '--cflags', package], env=env))
        # --static includes private dependencies and also works for shared PETSc.
        libs = shlex.split(output([pkgconfig, '--libs', '--static', package], env=env))
        version = output([pkgconfig, '--modversion', package], env=env)
        discovered = 'pkg-config:' + package
        prefix_text = output([pkgconfig, '--variable=prefix', package], env=env)
        if prefix_text:
            prefix = Path(prefix_text)
    elif root and (root / 'lib/petsc/conf/variables').is_file():
        # PETSc's own make variables handle in-place and prefix installations.
        with tempfile.TemporaryDirectory(prefix='slurry-petsc-flags-') as temp:
            temp = Path(temp)
            makefile = temp / 'flags.mk'
            makefile.write_text('include ' + str(root / 'lib/petsc/conf/variables') + '\n'
                                '.PHONY: all\nall:\n'
                                '\t$(file >' + str(temp / 'cflags') + ',$(PETSC_CC_INCLUDES))\n'
                                '\t$(file >' + str(temp / 'libs') + ',$(PETSC_LIB))\n'
                                '\t@true\n')
            output(['make', '--no-print-directory', '-s', '-f', str(makefile),
                    'PETSC_DIR=' + str(root), 'PETSC_ARCH=' + arch, 'all'], env=env)
            flags = shlex.split((temp / 'cflags').read_text())
            libs = shlex.split((temp / 'libs').read_text())
        version = 'header'
        discovered = 'PETSC_DIR/PETSC_ARCH'
    else:
        raise RuntimeError('PETSc not found. Submit build_slurry_cpu.sbatch, or set PETSC_DIR '
                           '(and PETSC_ARCH for an in-place build), or load an installed PETSc module.')
    if not libs:
        raise RuntimeError('PETSc discovery returned no link libraries')
    include_dirs = [Path(f[2:]) for f in flags if f.startswith('-I')]
    lib_dirs = [Path(f[2:]) for f in libs if f.startswith('-L')]
    dependencies = {}
    for directory in include_dirs:
        for name in ('petscconf.h', 'petscversion.h', 'petscsnes.h'):
            path = directory / name
            if path.is_file():
                dependencies[str(path.resolve())] = sha256(path)
    for directory in lib_dirs:
        # Record the actual library content so an in-place PETSc rebuild cannot
        # silently reuse an executable from the old build cache.
        found = False
        for name in ('libpetsc.so', 'libpetsc.dylib', 'libpetsc.a'):
            path = directory / name
            if path.is_file():
                dependencies[str(path.resolve())] = sha256(path)
                found = True
                break
        if found:
            break
    for item in libs:
        path = Path(item)
        if path.is_absolute() and path.is_file() and 'petsc' in path.name:
            dependencies[str(path.resolve())] = sha256(path)
    for directory in lib_dirs:
        if directory.is_dir():
            flag = '-Wl,-rpath,' + str(directory.resolve())
            if flag not in libs:
                libs.append(flag)
    return {'backend': 'petsc', 'discovery': discovered, 'version': version,
            'cflags': flags, 'libs': libs, 'dependencies_sha256': dependencies}


PROBE = r'''
#include <petscsnes.h>
#include <cstdio>
#if defined(PETSC_USE_COMPLEX) || !defined(PETSC_USE_REAL_DOUBLE)
#error Slurry requires a real double-precision PETSc build
#endif
#if PETSC_VERSION_LT(3,18,0)
#error Slurry requires PETSc 3.18 or later
#endif
#if defined(SLURRY_REQUIRE_MPI) && defined(PETSC_HAVE_MPIUNI)
#error The MPI slurry build requires PETSc built against the same MPI implementation
#endif
int main() {
  PetscErrorCode err = PetscInitializeNoArguments();
  if (err) return (int)err;
  SNES snes = nullptr;
  err = SNESCreate(PETSC_COMM_SELF, &snes);
  if (!err) err = SNESSetType(snes, SNESNEWTONLS);
  if (!err) err = SNESDestroy(&snes);
  char version[1024] = {0};
  PetscGetVersion(version, sizeof version);
  std::printf("PETSC=%s\n", version);
  std::printf("PETSC_HEADER=%d.%d.%d\n", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR);
#if !defined(PETSC_HAVE_MPIUNI)
  char mpi[MPI_MAX_LIBRARY_VERSION_STRING] = {0}; int length = 0;
  MPI_Get_library_version(mpi, &length);
  for (int i=0; i<length; ++i) if (mpi[i]=='\n' || mpi[i]=='\r') mpi[i]=' ';
  std::printf("MPI=%s\n", mpi);
#else
  std::printf("MPI=PETSc MPIUNI (serial)\n");
#endif
  PetscErrorCode finalError = PetscFinalize();
  return (int)(err ? err : finalError);
}
'''


def preflight(config, compiler, require_mpi=True):
    with tempfile.TemporaryDirectory(prefix='slurry-petsc-check-') as temp:
        source = Path(temp) / 'petsc_check.cpp'
        binary = Path(temp) / 'petsc_check'
        source.write_text(PROBE)
        command = [compiler, '-std=c++20'] + config['cflags']
        if require_mpi:
            command.append('-DSLURRY_REQUIRE_MPI')
        command += [str(source), '-o', str(binary)] + config['libs']
        try:
            output(command)
        except subprocess.CalledProcessError as error:
            raise RuntimeError('PETSc compile/link preflight failed. Use the same compiler '
                               'and MPI modules for PETSc, build, and run.\n' + error.output)
        launcher = mpi_launcher_prefix(compiler) if require_mpi else []
        try:
            report = output(launcher + [str(binary)])
        except subprocess.CalledProcessError as error:
            raise RuntimeError('PETSc compiled and linked, but its runtime check failed.\n'
                               'Launcher: ' + (' '.join(shlex.quote(x) for x in launcher) or 'serial direct execution') +
                               '\n' + error.output)
    config = dict(config)
    config['compiler'] = str(Path(compiler).absolute())
    config['compiler_version'] = output([compiler, '--version'])
    config['compiler_wrapper'] = wrapper_command(compiler)
    config['runtime_preflight'] = report
    # Keep the random temporary probe path out of the stable build fingerprint.
    config['runtime_launcher'] = launcher
    match = re.search(r'^PETSC_HEADER=(.+)$', report, re.MULTILINE)
    if match:
        config['version'] = match.group(1)
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default=os.environ.get('SLURRY_CXX', 'mpicxx'))
    parser.add_argument('--serial', action='store_true', help='Allow MPIUNI PETSc')
    parser.add_argument('--format', choices=('json', 'cflags', 'libs'), default='json')
    parser.add_argument('--no-probe', action='store_true', help='Discover only; for shell compile commands')
    args = parser.parse_args()
    config = discover()
    if not args.no_probe:
        compiler = shutil.which(args.compiler)
        if not compiler:
            raise RuntimeError('Compiler not found: ' + args.compiler)
        config = preflight(config, compiler, not args.serial)
    if args.format == 'json':
        print(json.dumps(config, indent=2))
    else:
        print(' '.join(shlex.quote(item) for item in config[args.format]))


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error))
