#!/usr/bin/env python3
import sys
from pathlib import Path

LIST_PATH = Path('analysis/single_enzyme_ligand_jobs.list')
ROOT = Path('jobs/human_test_set')

def main():
    if not LIST_PATH.exists():
        print(f"ERROR: {LIST_PATH} not found")
        sys.exit(1)
    removed = 0
    for name in [x.strip() for x in LIST_PATH.read_text().splitlines() if x.strip()]:
        mdir = ROOT / name / 'output_msa'
        if not mdir.is_dir():
            continue
        for p in mdir.glob('*.json'):
            try:
                p.unlink()
                removed += 1
            except Exception as e:
                print(f"WARN: failed to remove {p}: {e}")
    print(f"Removed {removed} JSON files under output_msa/ for listed jobs.")

if __name__ == '__main__':
    main()

