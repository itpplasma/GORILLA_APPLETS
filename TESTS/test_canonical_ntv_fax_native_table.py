"""Read an existing physical table with the native consumer at native samples."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

reader, table, executable, output = map(Path,sys.argv[1:5])
exports = list(map(Path,sys.argv[5:]))
spec = importlib.util.spec_from_file_location('neo2_fax',reader)
fax = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fax)
metadata = json.loads(table.with_suffix(table.suffix+'.json').read_text())
# Stream the large table hash without another table-sized allocation.
h = hashlib.sha256()
with table.open('rb') as f:
    for chunk in iter(lambda:f.read(1024*1024),b''):
        h.update(chunk)
if h.hexdigest() != metadata['sha256']:
    raise SystemExit('table checksum differs from its normalization receipt')
results = [fax.consumer_oracle(fax.read_fax(f),table,executable,
                              np.asarray(metadata['source_amplitudes'])) for f in exports]
output.write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
