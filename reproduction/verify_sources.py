from pathlib import Path
import csv
import hashlib

root = Path(__file__).resolve().parent
with (root / 'source_manifest.csv').open(newline='') as stream:
    rows = list(csv.DictReader(stream))
for row in rows:
    relative = Path(row['path'])
    if relative.is_absolute() or '..' in relative.parts:
        raise SystemExit('Unsafe manifest path')
    path = root / relative
    if not path.is_file() or path.stat().st_size != int(row['bytes']):
        raise SystemExit('Missing file or size mismatch: ' + row['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
        raise SystemExit('Checksum mismatch: ' + row['path'])
print('VERIFIED software-only source bundle:', len(rows), 'files; no fitting')
