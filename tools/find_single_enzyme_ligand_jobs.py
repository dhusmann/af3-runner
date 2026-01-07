#!/usr/bin/env python3
"""
Scan jobs/human_test_set for single-protein + ligand jobs whose directory
names match ENZYME-SAM or ENZYME-SAH, and verify their alphafold_input.json
contains exactly one protein and one ligand entry.

Outputs:
- TSV report to stdout
- Writes clean list to analysis/single_enzyme_ligand_jobs.list
- Also writes the TSV to analysis/single_enzyme_ligand_jobs.tsv
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path('jobs/human_test_set')
OUT_DIR = Path('analysis')
OUT_DIR.mkdir(parents=True, exist_ok=True)

name_re = re.compile(r'^[^-]+-(SAM|SAH)$', re.IGNORECASE)

if not ROOT.is_dir():
    print('ERROR: jobs/human_test_set not found', file=sys.stderr)
    sys.exit(1)

rows = []
for job_dir in sorted(p for p in ROOT.iterdir() if p.is_dir() and name_re.match(p.name)):
    jf = job_dir / 'alphafold_input.json'
    nprot = nlig = nother = 0
    note = ''
    if jf.exists():
        try:
            data = json.load(jf.open())
            for ent in data.get('sequences', []):
                if 'protein' in ent:
                    nprot += 1
                elif 'ligand' in ent:
                    nlig += 1
                else:
                    nother += 1
        except Exception as e:
            note = f'JSON_ERROR:{e.__class__.__name__}'
    status = (
        'OK_SINGLE'
        if jf.exists() and nprot == 1 and nlig == 1 and nother == 0
        else ('MISSING_JSON' if not jf.exists() else f'NOT_SINGLE(p{nprot},l{nlig},o{nother})')
    )
    rows.append((job_dir.name, status, nprot, nlig, nother, note))

# Emit TSV report
tsv_path = OUT_DIR / 'single_enzyme_ligand_jobs.tsv'
with tsv_path.open('w') as fh:
    fh.write('# job_name\tstatus\tn_proteins\tn_ligands\tn_other\tnotes\n')
    for r in rows:
        fh.write('\t'.join(map(str, r)) + '\n')

print('# job_name\tstatus\tn_proteins\tn_ligands\tn_other\tnotes')
for r in rows:
    print('\t'.join(map(str, r)))

# Write clean list
ok = [name for name, status, *_ in rows if status == 'OK_SINGLE']
lst_path = OUT_DIR / 'single_enzyme_ligand_jobs.list'
lst_path.write_text('\n'.join(ok) + ('\n' if ok else ''))

print('\n# Single-protein-ligand jobs (verified):')
for name in ok:
    print(name)
print(f'\nTotal OK single-protein-ligand: {len(ok)}')

