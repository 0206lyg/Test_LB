"""Archive only superseded project manuals during the existing common build."""
from datetime import datetime, timezone
import errno
import hashlib
import os
from pathlib import Path
import stat
import zipfile


OBSOLETE_DOCS = (
    'README_PETSC_BUILD.md',
    'README_PETSC_UPDATE_KO.md',
    'VALIDATION.md',
    'docs/graphite-local-adhesion.md',
    'docs/graphite-restart.md',
    'slurry/fitting/README_KO.md',
    'tests/petsc_contact/fixtures/README.md',
    'README_pure_gr_restart.txt',
)
_DIR_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW


def _directory(root_fd, parts, create=False):
    """Walk relative directory descriptors; never traverse symbolic links."""
    fd = os.dup(root_fd)
    try:
        for part in parts:
            if part in ('', '.', '..'):
                raise ValueError('Invalid documentation archive path component')
            if create:
                try:
                    os.mkdir(part, dir_fd=fd)
                except FileExistsError:
                    pass
            next_fd = os.open(part, _DIR_FLAGS, dir_fd=fd)
            os.close(fd)
            fd = next_fd
        return fd
    except BaseException:
        os.close(fd)
        raise


def _read_regular(parent_fd, name):
    fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=parent_fd)
    with os.fdopen(fd, 'rb') as stream:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode):
            raise ValueError('Obsolete documentation path is not a regular file: ' + name)
        return stream.read(), (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)


def _verify_archive(stream, records):
    stream.seek(0)
    with zipfile.ZipFile(stream, 'r') as archive:
        if archive.namelist() != [record['name'] for record in records]:
            raise RuntimeError('Documentation backup member verification failed')
        for record in records:
            if archive.read(record['name']) != record['data']:
                raise RuntimeError('Documentation backup content verification failed: ' + record['name'])


def archive_obsolete_docs(root):
    """Return archive path, removed names and skipped links; preserve all originals on archive failure."""
    root = Path(root)
    root_fd = os.open(str(root), _DIR_FLAGS)
    records, skipped = [], []
    backup_fd = None
    try:
        for name in OBSOLETE_DOCS:
            parts = name.split('/')
            parent_fd = None
            try:
                parent_fd = _directory(root_fd, parts[:-1])
                data, identity = _read_regular(parent_fd, parts[-1])
            except FileNotFoundError:
                if parent_fd is not None:
                    os.close(parent_fd)
                continue
            except (OSError, ValueError) as error:
                if parent_fd is not None:
                    os.close(parent_fd)
                if isinstance(error, OSError) and error.errno not in (errno.ELOOP, errno.ENOTDIR):
                    raise
                skipped.append(name)
                continue
            records.append({'name': name, 'data': data, 'identity': identity, 'parent_fd': parent_fd})
        if not records:
            return {'archive': None, 'removed': [], 'skipped': skipped}
        backup_fd = _directory(root_fd, ['build', 'documentation_backup'], create=True)
        content_hash = hashlib.sha256()
        for record in records:
            content_hash.update(record['name'].encode('utf-8') + b'\0')
            content_hash.update(hashlib.sha256(record['data']).digest())
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive_name = 'obsolete_documents_' + stamp + '_' + content_hash.hexdigest()[:12] + '.zip'
        partial_name = archive_name + '.partial'
        fd = os.open(partial_name, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                     0o600, dir_fd=backup_fd)
        try:
            with os.fdopen(fd, 'w+b') as stream:
                with zipfile.ZipFile(stream, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
                    for record in records:
                        archive.writestr(record['name'], record['data'])
                stream.flush()
                os.fsync(stream.fileno())
                _verify_archive(stream, records)
            os.rename(partial_name, archive_name, src_dir_fd=backup_fd, dst_dir_fd=backup_fd)
            os.fsync(backup_fd)
        except BaseException:
            try:
                os.unlink(partial_name, dir_fd=backup_fd)
            except FileNotFoundError:
                pass
            raise
        # Validate every source again before removing any. A concurrently edited
        # file stays in place; the verified archive already preserves the snapshot.
        for record in records:
            data, identity = _read_regular(record['parent_fd'], record['name'].split('/')[-1])
            if identity != record['identity'] or data != record['data']:
                raise RuntimeError('Documentation changed while being archived; originals retained: ' + record['name'])
        for record in records:
            os.unlink(record['name'].split('/')[-1], dir_fd=record['parent_fd'])
        return {'archive': str(root / 'build/documentation_backup' / archive_name),
                'removed': [record['name'] for record in records], 'skipped': skipped}
    finally:
        for record in records:
            os.close(record['parent_fd'])
        if backup_fd is not None:
            os.close(backup_fd)
        os.close(root_fd)
