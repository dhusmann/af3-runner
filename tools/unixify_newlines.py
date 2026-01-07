#!/usr/bin/env python3
"""Convert CSVs to Unix LF newlines in-place.
Targets: folding_jobs.csv, waiting_for_msa.csv, msa_array_jobs.csv
"""
from pathlib import Path

FILES = [
    Path('folding_jobs.csv'),
    Path('waiting_for_msa.csv'),
    Path('msa_array_jobs.csv'),
]

def unixify(p: Path):
    if not p.exists():
        return False
    b = p.read_bytes()
    nb = b.replace(b'\r\n', b'\n').replace(b'\r', b'\n')
    if nb != b:
        p.write_bytes(nb)
        return True
    return False

def main():
    changed = 0
    for f in FILES:
        if unixify(f):
            print(f"Fixed newlines: {f}")
            changed += 1
    print(f"Done. Files changed: {changed}")

if __name__ == '__main__':
    main()

