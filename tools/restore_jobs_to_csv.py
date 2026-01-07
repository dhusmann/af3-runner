#!/usr/bin/env python3
import csv
from pathlib import Path

CSV = Path('folding_jobs.csv')
LIST = Path('analysis/single_enzyme_ligand_jobs.list')

def read_csv_jobs(path: Path):
    if not path.exists():
        return []
    rows = []
    with path.open() as fh:
        r = csv.DictReader(fh)
        col = 'input_folder_name' if 'input_folder_name' in r.fieldnames else r.fieldnames[0]
        for row in r:
            val = row.get(col, '').strip()
            if val:
                rows.append(val)
    return rows

def write_csv_jobs(path: Path, jobs):
    # Force Unix line endings to avoid CRLF issues downstream
    with path.open('w', newline='') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['input_folder_name'])
        for j in jobs:
            w.writerow([j])

def main():
    if not LIST.exists():
        print('ERROR: analysis/single_enzyme_ligand_jobs.list not found')
        return
    current = set(read_csv_jobs(CSV))
    to_add = [j.strip() for j in LIST.read_text().splitlines() if j.strip()]
    merged = list(sorted(current.union(to_add))) if current else to_add
    write_csv_jobs(CSV, merged)
    print(f'folding_jobs.csv updated. Total entries: {len(merged)}')

if __name__ == '__main__':
    main()
