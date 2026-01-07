#!/usr/bin/env python3
import sys
import shutil
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        print("Usage: remove_outputs_from_list.py <path-to-job-list>")
        sys.exit(1)
    list_path = Path(sys.argv[1])
    base = Path('jobs/human_test_set')
    if not list_path.exists():
        print(f"ERROR: list file not found: {list_path}")
        sys.exit(1)
    removed = 0
    skipped_missing = 0
    total = 0
    lines = [ln.strip().strip('"').strip("'\t ") for ln in list_path.read_text().splitlines()]
    for name in lines:
        if not name or name.startswith('#'):
            continue
        total += 1
        outdir = base / name / 'output'
        if outdir.is_dir():
            try:
                shutil.rmtree(outdir)
                removed += 1
            except Exception as e:
                print(f"WARN: failed to remove {outdir}: {e}")
        else:
            skipped_missing += 1
    print(f"Processed {total} job names from {list_path}")
    print(f"Removed output/ dirs: {removed}")
    print(f"No output/ present:   {skipped_missing}")

if __name__ == '__main__':
    main()

