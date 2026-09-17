#!/usr/bin/env python3
"""Recover the exact original bytes after removing integration-only edits."""
import hashlib
import json
from pathlib import Path
import sys

def main():
    base = Path(__file__).resolve().parents[2]
    manifest = json.loads((base / 'slurry/source_manifest.json').read_text())
    failures = []
    for item in manifest['files']:
        data = (base / item['destination']).read_bytes()
        if item['namespace']:
            text = data.decode('utf-8')
            for marker in (item['opening'], item['closing']):
                if text.count(marker) != 1:
                    failures.append(item['destination'] + ': integration marker changed'); break
                text = text.replace(marker, '', 1)
            if item['entry']:
                if text.count('int runCase(') != 1: failures.append(item['destination'] + ': entry point changed')
                text = text.replace('int runCase(', 'int main(', 1)
            for before, after in item.get('guards', []):
                text = text.replace(after, before)
            data = text.encode('utf-8')
        if hashlib.sha256(data).hexdigest() != item['source_sha256']:
            failures.append(item['destination'] + ': differs from pinned original')
    if failures:
        print('\n'.join(failures)); return 1
    print('PASS: all %d original source/data files preserved; only C++ namespace/entry-point/header-guard edits.' % len(manifest['files']))
    print('Source commit: ' + manifest['commit'])
    return 0

if __name__ == '__main__': sys.exit(main())
