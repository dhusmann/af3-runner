#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

def read_csv_jobs(path: Path):
    if not path.exists():
        return []
    rows = []
    with path.open(newline='') as fh:
        r = csv.DictReader(fh)
        # be tolerant of header name
        col = 'input_folder_name' if 'input_folder_name' in r.fieldnames else r.fieldnames[0]
        for row in r:
            val = (row.get(col) or '').strip().strip('"')
            if val:
                rows.append(val)
    return rows

def write_csv_jobs(path: Path, jobs):
    # Always write LF to avoid CRLF issues downstream
    with path.open('w', newline='') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['input_folder_name'])
        for j in jobs:
            w.writerow([j])

def read_list(list_path: Path):
    names = []
    for ln in list_path.read_text().splitlines():
        ln = ln.strip().strip('"').strip("'\t ")
        if not ln or ln.startswith('#'):
            continue
        names.append(ln)
    return names

def main():
    if len(sys.argv) < 2:
        print('Usage: restore_jobs_from_list.py <path-to-job-list>')
        sys.exit(1)
    list_path = Path(sys.argv[1])
    csv_path = Path('folding_jobs.csv')
    if not list_path.exists():
        print(f'ERROR: list file not found: {list_path}')
        sys.exit(1)
    current = set(read_csv_jobs(csv_path))
    incoming = set(read_list(list_path))
    merged = sorted(current.union(incoming))
    write_csv_jobs(csv_path, merged)
    print(f'folding_jobs.csv updated. Total entries: {len(merged)}')
    print(f'Added from list (including already present): {len(incoming)}')

if __name__ == '__main__':
    main()

