"""The single build archives only explicitly retired manuals, without data loss."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
sys.path.insert(0, str(TOOLS))
import build_slurry
import consolidate_docs as migration


class DocumentationMigrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='slurry-doc-migration-')
        self.root = Path(self.temporary.name) / 'OpenLB'
        self.root.mkdir()

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, name, data):
        target = self.root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        return target

    def test_exact_archive_idempotence_and_unrelated_files_preserved(self):
        expected = {}
        for index, name in enumerate(migration.OBSOLETE_DOCS):
            data = ('원본 %d\r\n' % index).encode('utf-8') + b'\x00\xff'
            self.write(name, data)
            expected[name] = data
        preserve = ['README.md', 'slurry/cases/pure_gr.json', 'olb-1.9r0/README.md',
                    'olb-1.9r0/LICENSE', 'vendor/README.md', 'docs/my_notes.md']
        for name in preserve:
            self.write(name, b'keep exactly\r\n')
        result = migration.archive_obsolete_docs(self.root)
        self.assertEqual(result['removed'], list(migration.OBSOLETE_DOCS))
        self.assertEqual(result['skipped'], [])
        with zipfile.ZipFile(result['archive']) as archive:
            self.assertEqual(archive.namelist(), list(migration.OBSOLETE_DOCS))
            self.assertEqual({name: archive.read(name) for name in archive.namelist()}, expected)
        for name in migration.OBSOLETE_DOCS:
            self.assertFalse((self.root / name).exists())
        for name in preserve:
            self.assertEqual((self.root / name).read_bytes(), b'keep exactly\r\n')
        self.assertIsNone(migration.archive_obsolete_docs(self.root)['archive'])
        self.assertEqual(len(list((self.root / 'build/documentation_backup').glob('*.zip'))), 1)

    def test_verification_failure_leaves_every_original(self):
        files = {name: self.write(name, name.encode()) for name in migration.OBSOLETE_DOCS}
        with patch.object(migration, '_verify_archive', side_effect=RuntimeError('verification failed')):
            with self.assertRaisesRegex(RuntimeError, 'verification failed'):
                migration.archive_obsolete_docs(self.root)
        for name, path in files.items():
            self.assertEqual(path.read_bytes(), name.encode())
        self.assertEqual(list((self.root / 'build/documentation_backup').iterdir()), [])

    def test_changed_source_preserves_all_originals_and_verified_snapshot(self):
        first = self.write('README_PETSC_BUILD.md', b'original build notes')
        second = self.write('VALIDATION.md', b'original validation')
        verify = migration._verify_archive
        def change_after_verify(stream, records):
            verify(stream, records)
            second.write_bytes(b'new user note')
        with patch.object(migration, '_verify_archive', side_effect=change_after_verify):
            with self.assertRaisesRegex(RuntimeError, 'Documentation changed'):
                migration.archive_obsolete_docs(self.root)
        self.assertEqual(first.read_bytes(), b'original build notes')
        self.assertEqual(second.read_bytes(), b'new user note')
        archives = list((self.root / 'build/documentation_backup').glob('*.zip'))
        self.assertEqual(len(archives), 1)
        with zipfile.ZipFile(archives[0]) as archive:
            self.assertEqual(archive.read('VALIDATION.md'), b'original validation')

    def test_file_and_parent_symlinks_are_not_followed(self):
        outside = Path(self.temporary.name) / 'outside'
        outside.mkdir()
        target = outside / 'graphite-local-adhesion.md'
        target.write_bytes(b'outside content')
        (self.root / 'docs').symlink_to(outside, target_is_directory=True)
        (self.root / 'VALIDATION.md').symlink_to(target)
        self.write('README_PETSC_BUILD.md', b'normal document')
        result = migration.archive_obsolete_docs(self.root)
        self.assertEqual(result['removed'], ['README_PETSC_BUILD.md'])
        self.assertIn('VALIDATION.md', result['skipped'])
        self.assertIn('docs/graphite-local-adhesion.md', result['skipped'])
        self.assertTrue((self.root / 'docs').is_symlink())
        self.assertTrue((self.root / 'VALIDATION.md').is_symlink())
        self.assertEqual(target.read_bytes(), b'outside content')

    def test_backup_parent_symlink_rejects_before_removing_sources(self):
        source = self.write('VALIDATION.md', b'original')
        outside = Path(self.temporary.name) / 'outside'
        outside.mkdir()
        (self.root / 'build').symlink_to(outside, target_is_directory=True)
        with self.assertRaises(OSError):
            migration.archive_obsolete_docs(self.root)
        self.assertEqual(source.read_bytes(), b'original')
        self.assertEqual(list(outside.iterdir()), [])

    def test_existing_build_entry_migrates_even_on_cache_hit(self):
        self.write('README_pure_gr_restart.txt', b'previous package instructions')
        olb = self.root / 'olb-1.9r0'
        self.write('olb-1.9r0/src/case/case.h', b'header')
        self.write('build/slurry/releases/cache/slurry', b'executable')
        self.write('build/slurry/releases/cache/build_manifest.json',
                   json.dumps({'executable_sha256': 'hash'}).encode())
        with patch.object(build_slurry, 'BASE', self.root), patch.object(build_slurry, 'OLB', olb), \
                patch.object(build_slurry.shutil, 'which', return_value='/compiler/g++') as compiler, \
                patch.object(build_slurry, 'tool_output', return_value='compiler version'), \
                patch.object(build_slurry, 'digest', return_value='hash'), \
                patch.object(build_slurry, 'source_inputs', return_value={}), \
                patch.object(build_slurry, 'fingerprint', return_value='cache'), \
                patch.object(build_slurry, 'publish') as publish, \
                patch.object(build_slurry.subprocess, 'run') as compile_process, \
                patch.object(sys, 'argv', ['build_slurry.py', '--mode', 'serial', '--particle-solver', 'legacy']):
            self.assertEqual(build_slurry.main(), 0)
        compiler.assert_called_once_with('g++')
        publish.assert_called_once()
        compile_process.assert_not_called()
        self.assertFalse((self.root / 'README_pure_gr_restart.txt').exists())
        archives = list((self.root / 'build/documentation_backup').glob('*.zip'))
        self.assertEqual(len(archives), 1)
        with zipfile.ZipFile(archives[0]) as archive:
            self.assertEqual(archive.read('README_pure_gr_restart.txt'), b'previous package instructions')


if __name__ == '__main__':
    unittest.main()
