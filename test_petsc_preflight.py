"""Regressions for one-rank MPI launch and reuse after a failed install check."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
sys.path.insert(0, str(TOOLS))
import petsc_config
import setup_petsc


class PreflightTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='slurry-launch-test-')
        self.root = Path(self.temporary.name)
        self.compiler = self.root / 'mpi/bin/mpicxx'
        self.compiler.parent.mkdir(parents=True)
        self.launcher = self.compiler.with_name('mpirun')
        for executable in (self.compiler, self.launcher):
            executable.write_text('#!/bin/sh\nexit 0\n')
            executable.chmod(0o755)
        self.config = {'cflags': [], 'libs': [], 'version': 'header'}

    def tearDown(self):
        self.temporary.cleanup()

    def fake_output(self, argv):
        if argv[0] == str(self.compiler):
            if argv[-1] == '--version':
                return 'GCC test'
            if argv[-1] in ('--showme', '-show'):
                return 'g++ -lmpi'
            self.assertIn('-DSLURRY_REQUIRE_MPI', argv)
            return ''
        if argv == [str(self.launcher), '-np', '1', argv[-1]]:
            self.assertEqual(Path(argv[-1]).name, 'petsc_check')
            return 'PETSC_HEADER=3.25.5\nMPI=test OpenMPI'
        raise subprocess.CalledProcessError(1, argv, 'MPI_Init_thread: direct launch fails under Slurm')

    def test_compiler_sibling_wins_over_unrelated_path_launcher(self):
        with patch.object(petsc_config.shutil, 'which', return_value='/wrong/mpi/bin/mpirun'):
            self.assertEqual(petsc_config.mpi_launcher_prefix(str(self.compiler)),
                             [str(self.launcher), '-np', '1'])

    def test_mpiexec_sibling_fallback(self):
        self.launcher.rename(self.launcher.with_name('mpiexec'))
        self.assertEqual(petsc_config.mpi_launcher_prefix(str(self.compiler))[0],
                         str(self.launcher.with_name('mpiexec')))

    def test_missing_launcher_is_actionable(self):
        self.launcher.unlink()
        with patch.object(petsc_config.shutil, 'which', return_value=None):
            with self.assertRaisesRegex(RuntimeError, 'MPI launcher not found'):
                petsc_config.mpi_launcher_prefix(str(self.compiler))

    def test_slurm_probe_uses_launcher_and_fingerprint_is_stable(self):
        with patch.object(petsc_config, 'output', side_effect=self.fake_output), \
                patch.dict(os.environ, {'SLURM_JOB_ID': '22944246', 'SLURM_NTASKS': '32'}):
            first = petsc_config.preflight(self.config, str(self.compiler))
            second = petsc_config.preflight(self.config, str(self.compiler))
        self.assertEqual(first, second)
        self.assertEqual(first['version'], '3.25.5')
        self.assertEqual(first['runtime_launcher'], [str(self.launcher), '-np', '1'])

    def test_runtime_failure_is_not_reported_as_compile_failure(self):
        def failed_runtime(argv):
            if argv[0] == str(self.launcher):
                raise subprocess.CalledProcessError(1, argv, 'MPI_Init_thread failed')
            return self.fake_output(argv)
        with patch.object(petsc_config, 'output', side_effect=failed_runtime):
            with self.assertRaisesRegex(RuntimeError, 'compiled and linked, but its runtime check failed'):
                petsc_config.preflight(self.config, str(self.compiler))

    def test_compile_failure_does_not_launch(self):
        with patch.object(petsc_config, 'output', side_effect=subprocess.CalledProcessError(1, ['compiler'], 'link failed')), \
                patch.object(petsc_config, 'mpi_launcher_prefix') as launcher:
            with self.assertRaisesRegex(RuntimeError, 'compile/link preflight failed'):
                petsc_config.preflight(self.config, str(self.compiler))
            launcher.assert_not_called()

    def test_serial_probe_does_not_require_mpi_launcher(self):
        with patch.object(petsc_config, 'output', return_value='PETSC_HEADER=3.25.5\nMPI=MPIUNI'), \
                patch.object(petsc_config, 'mpi_launcher_prefix') as launcher:
            result = petsc_config.preflight(self.config, str(self.compiler), False)
        self.assertEqual(result['runtime_launcher'], [])
        launcher.assert_not_called()

    def test_completed_install_is_reused_without_downloading_or_make(self):
        ccompiler = self.compiler.with_name('mpicc')
        identity = {'cxx': str(self.compiler), 'cc': str(ccompiler),
                    'wrapper': 'g++ -lmpi', 'version': 'GCC test'}
        fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:12]
        prefix = self.root / 'installs' / setup_petsc.VERSION / fingerprint / 'install'
        (prefix / 'include').mkdir(parents=True)
        (prefix / 'include/petscconf.h').write_text('/* already installed */\n')
        config = {'version': '3.25.5', 'runtime_launcher': [str(self.launcher), '-np', '1']}
        def find(name):
            return str(ccompiler if name == 'mpicc' else self.compiler)
        with patch.dict(os.environ, {'SLURRY_CXX': 'mpicxx', 'SLURRY_CC': 'mpicc'}, clear=True), \
                patch.object(setup_petsc, 'BASE', self.root / 'OpenLB'), \
                patch.object(setup_petsc.shutil, 'which', side_effect=find), \
                patch.object(setup_petsc, 'wrapper_command', return_value='g++ -lmpi'), \
                patch.object(setup_petsc, 'output', return_value='GCC test'), \
                patch.object(setup_petsc, 'check_existing', side_effect=[RuntimeError('PETSc not found.'), config]), \
                patch.object(setup_petsc.urllib.request, 'urlretrieve') as download, \
                patch.object(setup_petsc.subprocess, 'run') as build:
            setup_petsc.main(['--prefix-root', str(self.root / 'installs')])
            self.assertEqual(os.environ['PETSC_DIR'], str(prefix))
            download.assert_not_called()
            build.assert_not_called()
        self.assertTrue((self.root / 'OpenLB/build/petsc/env.sh').is_file())


if __name__ == '__main__':
    unittest.main()
