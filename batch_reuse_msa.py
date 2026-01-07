#!/usr/bin/env python3
"""
batch_reuse_msa.py  (v7.2, 2025-08-04)
────────────────────────────────────
Fast MSA re-use *plus* triage of missing jobs into:

  • msa_array_jobs.csv       ← first MISS for each entry-pair
  • waiting_for_msa.csv      ← later MISSes that depend on the above

Run:
  python batch_reuse_msa.py --dry-run   # preview (no writes, no copies)
  python batch_reuse_msa.py             # copy MSAs & rewrite the CSVs
"""
import argparse, csv, glob, json, pathlib, re, sys

# ─── settings ────────────────────────────────────────────────────────────
JOBS_DIR  = pathlib.Path("jobs")
CSV_MAIN  = pathlib.Path("folding_jobs.csv")  # Fixed: was hts_folding_jobs.csv
CSV_MSA   = pathlib.Path("msa_array_jobs.csv")
CSV_WAIT  = pathlib.Path("waiting_for_msa.csv")
CORE_KEYS = {"id","sequence","modifications"}

parser = argparse.ArgumentParser()
parser.add_argument("-n","--dry-run",action="store_true",
                    help="preview actions without writing files")
DRY = parser.parse_args().dry_run
# ─── helper functions ────────────────────────────────────────────────────
def parse_entries(job:str):
    e1, rest = job.split("-",1)
    if "_noM" in rest:
        e2 = re.match(r"([^_]*?_noM)",rest).group(1)
    else:
        e2 = re.match(r"([^_-]+)",rest).group(1)
    return e1, e2

def base_key(entry1:str, entry2:str) -> str:
    """strip PTM underscore (_K.., _H.., _Q..) from entry2"""
    base2 = re.sub(r'_[KHQ]\d.*','',entry2)   # keep ..._noM
    return f"{entry1}-{base2}"

def chain_map(seq):
    return {e["protein"]["id"]:e["protein"] for e in seq if "protein" in e}

def collapse(txt:str)->str:
    pat=r'("templateIndices"\s*:\s*\[)([\s0-9,]+?)(\])'
    return re.sub(pat,
                  lambda m:f'{m.group(1)}{m.group(2).replace("\\n","").replace(" ","")}{m.group(3)}',
                  txt, flags=re.DOTALL)

def merge_json(msa,fresh)->str:
    s=json.loads(msa.read_text()); d=json.loads(fresh.read_text())
    sc,dc = chain_map(s["sequences"]), chain_map(d["sequences"])
    miss=[cid for cid in dc if cid not in sc]
    if miss: raise ValueError(f"chains {miss} missing in {msa.name}")
    for ent in d["sequences"]:
        if "protein" not in ent: continue
        cid=ent["protein"]["id"]
        for k,v in sc[cid].items():
            if k not in CORE_KEYS: ent["protein"][k]=v
    d["version"]=max(d.get("version",1),s.get("version",1))
    return collapse(json.dumps(d,indent=2,separators=(",",": ")))

# ─── restrict search roots (top + one level, no "-") ────────────────────
TOP = [JOBS_DIR]
SUB = [p for p in JOBS_DIR.iterdir() if p.is_dir() and "-" not in p.name]
ROOTS = TOP + SUB
# ─── locate current job dir ──────────────────────────────────────────────
def find_job_dir(job):
    for root in ROOTS:
        cand=root/job
        if cand.is_dir() and (cand/"alphafold_input.json").is_file():
            return cand
    return None
def is_single_protein_ligand(job_dir:pathlib.Path) -> bool:
    try:
        d=json.loads((job_dir/"alphafold_input.json").read_text())
    except Exception:
        return False
    nprot=sum(1 for e in d.get("sequences",[]) if "protein" in e)
    nlig=sum(1 for e in d.get("sequences",[]) if "ligand" in e)
    return nprot==1 and nlig==1
# ─── find MSA candidates ────────────────────────────────────────────────
def find_msa_candidates(e1,e2,exclude:pathlib.Path):
    c=[]
    for root in ROOTS:
        # Use a pattern that matches both direct files and files in subdirs
        # The ** pattern with recursive=True will match 0 or more directories
        base_path = root / f"{e1}-{e2}*" / "output_msa"
        
        # Find all *_data.json files under output_msa (at any depth)
        for job_dir in glob.glob(str(base_path).replace("/output_msa", "")):
            msa_dir = pathlib.Path(job_dir) / "output_msa"
            if msa_dir.exists():
                # Find all *.json files in this directory tree
                for p in msa_dir.rglob("*.json"):
                    if exclude and exclude in p.parents: 
                        continue
                    c.append(p)
    return sorted(c)

# ─── main pass ───────────────────────────────────────────────────────────
msa_array, waiting = [], []
seen_base, base_done=set(), set()

stats={"copy":0,"skip":0,"warn":0}

if not CSV_MAIN.is_file(): sys.exit(f"missing {CSV_MAIN}")

# Check if CSV has header
with CSV_MAIN.open() as fh:
    first_line = fh.readline().strip()
    if not first_line or (not "input_folder_name" in first_line and not "folder" in first_line):
        print(f"ERROR: {CSV_MAIN} is missing header or is malformed")
        print(f"First line: '{first_line}'")
        print("Expected header: 'input_folder_name' or 'folder'")
        sys.exit(1)

with CSV_MAIN.open() as fh:
    reader = csv.DictReader(fh)
    # Handle both possible header names
    if "input_folder_name" in reader.fieldnames:
        field_name = "input_folder_name"
    elif "folder" in reader.fieldnames:
        field_name = "folder"
    else:
        print(f"ERROR: CSV header must contain 'input_folder_name' or 'folder'")
        print(f"Found headers: {reader.fieldnames}")
        sys.exit(1)
    
    rows = list(reader)

for row in rows:
    job=row[field_name].strip()
    e1,e2=parse_entries(job)
    bkey=base_key(e1,e2)

    job_dir=find_job_dir(job)
    if not job_dir:
        print(f"[WARN] job dir missing   : {job}")
        stats["warn"]+=1; continue
    # For single-protein+ligand jobs, do NOT reuse/copy MSAs; let AF3 run pipeline
    if is_single_protein_ligand(job_dir):
        stats["skip"]+=1
        base_done.add(bkey)
        continue

    out_dir=job_dir/"output_msa"
    if out_dir.is_dir():
        stats["skip"]+=1
        base_done.add(bkey)
        continue  # already has MSA

    # try to reuse an existing MSA
    cands = find_msa_candidates(e1, e2, exclude=job_dir)
    fresh = job_dir / "alphafold_input.json"

    if cands:                                     # <── NEW unified branch
        if DRY:
            print(f"[COPY] {job:30} ← {cands[0].parents[2].name}")
            base_done.add(bkey)                   # mark as satisfied
            continue
        try:
            merged = merge_json(cands[0], fresh)
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "alphafold_input_with_msa.json").write_text(merged)
            print(f"[COPY] {job:30} ← {cands[0].parents[2].name}")
            stats["copy"] += 1
            base_done.add(bkey)
            continue
        except Exception as e:
            print(f"[ERR ] {job}: {e}")
            stats["warn"] += 1

    # still MISS — triage
    if bkey in base_done or bkey in seen_base:
        waiting.append(job)
    else:
        msa_array.append(job)
        seen_base.add(bkey)

# ─── write / preview the two auxiliary CSVs ──────────────────────────────
def write_csv(path, items):
    # Write with Unix line endings to avoid CRLF issues in shell readers
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["input_folder_name"]) 
        w.writerows([[x] for x in items])

if DRY:
    print("\nWould write msa_array_jobs.csv:")
    print("\n".join("  "+x for x in msa_array) or "  (none)")
    print("\nWould write waiting_for_msa.csv:")
    print("\n".join("  "+x for x in waiting) or "  (none)")
else:
    write_csv(CSV_MSA,  msa_array)
    write_csv(CSV_WAIT, waiting)
    print(f"\nWrote {CSV_MSA}  ({len(msa_array)} jobs)")
    print(f"Wrote {CSV_WAIT} ({len(waiting)} jobs)")

# ─── final summary ───────────────────────────────────────────────────────
print("\nSummary\n-------")
print(f"copied MSAs : {stats['copy']}")
print(f"skipped     : {stats['skip']} (already had MSA)")
print(f"warnings    : {stats['warn']}")
print(f"to_generate : {len(msa_array)}  → msa_array_jobs.csv")
print(f"waiting     : {len(waiting)}   → waiting_for_msa.csv")
