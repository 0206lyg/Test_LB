"""Build/run bookkeeping only; no physical models or unit conversion."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

BASE = Path(__file__).resolve().parents[2]
OLB = BASE / 'olb-1.9r0'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def write_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)

def source_inputs():
    files = [OLB / n for n in ('config.mk', 'rules.mk')]
    files += [BASE / 'slurry/Makefile', BASE / 'slurry/engines.json', BASE / 'slurry/tools/build_slurry.py']
    for directory in (OLB / 'src', OLB / 'examples/slurry'):
        files.extend(p for p in directory.rglob('*') if p.is_file() and p.suffix in ('.h', '.hh', '.hpp', '.cpp', '.c'))
    return {str(p.relative_to(BASE)): digest(p) for p in sorted(set(files))}

def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode('utf-8')).hexdigest()

def tool_output(argv):
    return subprocess.check_output(argv, stderr=subprocess.STDOUT, universal_newlines=True).strip()

def git_state():
    try:
        return {'commit': tool_output(['git', '-C', str(BASE), 'rev-parse', 'HEAD']),
                'status': tool_output(['git', '-C', str(BASE), 'status', '--porcelain'])}
    except (OSError, subprocess.CalledProcessError):
        return None

def build_record(executable):
    executable = Path(executable).resolve()
    manifest = executable.parent / 'build_manifest.json'
    if not manifest.is_file():
        raise ValueError('Missing build_manifest.json beside executable; run the common build first')
    record = json.loads(manifest.read_text(encoding='utf-8'))
    if digest(executable) != record['executable_sha256']:
        raise ValueError('Executable checksum does not match its completed build')
    if source_inputs() != record['source_inputs']:
        raise ValueError('C++ source/build configuration changed. Run build_slurry_cpu.sbatch again')
    return record
