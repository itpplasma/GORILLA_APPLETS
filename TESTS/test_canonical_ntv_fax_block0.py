"""Bitwise legacy block-0 compatibility with a supplied #80 native baseline."""
from pathlib import Path
import hashlib
import sys

baseline, current = map(Path, sys.argv[1:])
for case in ('fax_off', 'fax_zero', 'fax_on'):
    name = Path(case)/'canonical_ntv_markers.dat'
    a, b = (baseline/name).read_bytes(), (current/name).read_bytes()
    if a != b:
        raise SystemExit(f'FAIL: block-0 compatibility: {case}')
    print(f'{case}: bitwise identical, sha256={hashlib.sha256(b).hexdigest()}')
