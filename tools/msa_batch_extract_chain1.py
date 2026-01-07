#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, os, re, sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Defaults can be overridden on the CLI
ROOT_DEFAULT = Path("/scratch/groups/ogozani/alphafold3/jobs/human_test_set")
CSV_DEFAULT = Path("folding_jobs_nsd2i.csv")

PAIR_CONTAINER_KEYS = {"paired_msa", "paired_msas", "complex_msa", "interaction_msa"}
DROP_EXACT_KEYS = {"pairedMsa"} | PAIR_CONTAINER_KEYS

@dataclass
class JobContext:
    job_name: str
    dest_dir: Path
    dest_input_json: Path
    donor_dir: Path
    donor_with_msa_json: Path

def read_lines(path: Path) -> List[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]

def load_json(path: Path) -> Any:
    with path.open("r") as f:
        return json.load(f)

def safe_write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)  # atomic

def is_aa_seq(s: str) -> bool:
    return bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", s)) and len(s) >= 30

def find_first_sequence(obj: Any) -> Optional[str]:
    best = None
    def visit(x):
        nonlocal best
        if isinstance(x, dict):
            for v in x.values(): visit(v)
        elif isinstance(x, list):
            for v in x: visit(v)
        elif isinstance(x, str):
            if is_aa_seq(x) and (best is None or len(x) > len(best)): best = x
    visit(obj)
    return best

def find_ligand_value(obj: Any) -> Optional[Any]:
    # Prefer ligand under sequences if present
    if isinstance(obj, dict) and isinstance(obj.get("sequences"), list):
        for el in obj["sequences"]:
            if isinstance(el, dict) and "ligand" in el:
                return el["ligand"]
    # Fallback: first ligand/ligands anywhere
    q=[obj]
    while q:
        cur=q.pop(0)
        if isinstance(cur, dict):
            if "ligand" in cur: return cur.get("ligand")
            if "ligands" in cur: return cur.get("ligands")
            q.extend(cur.values())
        elif isinstance(cur, list):
            q.extend(cur)
    return None

def listdirs(p: Path) -> List[Path]:
    return [x for x in p.iterdir() if x.is_dir()]

def donors_by_prefix(root: Path) -> Dict[str, List[Path]]:
    byp: Dict[str, List[Path]] = {}
    for d in listdirs(root):
        name = d.name
        if "-" not in name:
            continue
        if (d / "output_msa" / "alphafold_input_with_msa.json").is_file():
            byp.setdefault(name.split("-")[0], []).append(d)
    for lst in byp.values():
        lst.sort(key=lambda p: p.name)
    return byp

def find_donor(prefix_index: Dict[str, List[Path]], prefix: str, exclude: str) -> Optional[Path]:
    c = [p for p in prefix_index.get(prefix, []) if p.name != exclude]
    return c[0] if c else None

def looks_chain_like_element(x: Any) -> bool:
    if not isinstance(x, dict): return False
    keys=set(x.keys())
    hints={"msa","templates","sequence","seq","chain_id","alignments","unpairedMsa","pairedMsa"}
    return bool(keys & hints)

def is_chain_like_list(lst: List[Any]) -> bool:
    if len(lst)<2: return False
    dicts=[x for x in lst if isinstance(x, dict)]
    if len(dicts) < max(2, len(lst)//2): return False
    return any(looks_chain_like_element(x) for x in dicts)

def strip_paired_blocks(d: Any) -> Any:
    if isinstance(d, dict):
        out={}
        for k,v in d.items():
            # Drop explicit paired* keys/containers; never drop unpairedMsa
            if k in DROP_EXACT_KEYS or k.lower().startswith("paired"):
                continue
            out[k]=strip_paired_blocks(v)
        return out
    if isinstance(d, list):
        return [strip_paired_blocks(x) for x in d]
    return d

def choose_chain_index_by_sequence(donor: Any, dest_seq: Optional[str]) -> int:
    if not dest_seq: return 0
    candidates=[]
    if isinstance(donor, dict):
        for k in ("chains","sequences","entities","monomers","inputs","targets"):
            v=donor.get(k)
            if isinstance(v,list) and is_chain_like_list(v): candidates.append(v)
    if not candidates:
        def scan(o):
            if isinstance(o, dict):
                for vv in o.values(): scan(vv)
            elif isinstance(o, list):
                if is_chain_like_list(o): candidates.append(o)
                else:
                    for vv in o: scan(vv)
        scan(donor)
    def seq_from_elem(e):
        if isinstance(e, dict):
            for key in ("sequence","seq","query_sequence","target_sequence"):
                s=e.get(key)
                if isinstance(s,str) and is_aa_seq(s): return s
        return None
    best_idx,best_sc=0,-1
    for lst in candidates:
        for idx in range(min(2,len(lst))):
            s=seq_from_elem(lst[idx]); sc=0
            if s:
                if dest_seq in s or s in dest_seq: sc=3
                else:
                    n=min(len(s),len(dest_seq))
                    if n:
                        m=sum(1 for a,b in zip(s[:n],dest_seq[:n]) if a==b)
                        sc=int(100*m/max(1,n))
            if sc>best_sc: best_idx,best_sc=idx,sc
    return best_idx

def prune_chain_dimension(obj: Any, keep_index: int, chain_keys=("A","B")) -> Any:
    if isinstance(obj, dict):
        keys=set(obj.keys())
        ab=keys.issuperset(set(chain_keys))
        one_two=keys.issuperset({"1","2"}) or keys.issuperset({"chain_1","chain_2"})
        if ab or one_two:
            key_keep=chain_keys[keep_index] if ab else ("1" if keep_index==0 else "2")
            if one_two and ("chain_1" in keys or "chain_2" in keys):
                key_keep="chain_1" if keep_index==0 else "chain_2"
            if key_keep in obj: return prune_chain_dimension(obj[key_keep], keep_index, chain_keys)
        out={}
        for k,v in obj.items():
            # Remove paired everywhere
            if k in DROP_EXACT_KEYS or k.lower().startswith("paired"):
                continue
            # Row-level dicts like {"seq":{"A": "...", "B":"..."}} → keep only selected chain
            if k=="seq" and isinstance(v, dict):
                for cand in (chain_keys[keep_index], str(keep_index+1)):
                    if cand in v:
                        out[k]=prune_chain_dimension(v[cand], keep_index, chain_keys); break
                else:
                    out[k]=prune_chain_dimension(v, keep_index, chain_keys)
                continue
            out[k]=prune_chain_dimension(v, keep_index, chain_keys)
        return out
    if isinstance(obj, list):
        # Special-case: sequences array of {protein:{...}} / {ligand:{...}}
        if obj and all(isinstance(x, dict) for x in obj) and any(("protein" in x or "ligand" in x) for x in obj):
            out_list=[]
            prot_seen=0
            for el in obj:
                if "ligand" in el:
                    out_list.append(prune_chain_dimension(el, keep_index, chain_keys))
                elif "protein" in el:
                    if prot_seen==keep_index:
                        out_list.append(prune_chain_dimension(el, keep_index, chain_keys))
                    prot_seen+=1
            return out_list
        # Generic chain-like list
        if is_chain_like_list(obj) and len(obj)>keep_index:
            return prune_chain_dimension(obj[keep_index], keep_index, chain_keys)
        return [prune_chain_dimension(x, keep_index, chain_keys) for x in obj]
    return obj

def override_ligand(out_json: Any, ligand_value: Any) -> Any:
    if ligand_value is None: return out_json
    def walk(x):
        if isinstance(x, dict):
            y={}
            for k,v in x.items():
                if k=="ligand": y[k]=ligand_value
                elif k=="ligands": y[k]=ligand_value if isinstance(ligand_value,list) else [ligand_value]
                else: y[k]=walk(v)
            return y
        if isinstance(x, list): return [walk(v) for v in x]
        return x
    return walk(out_json)

def output_is_complete(out_json: Any, dest_json: Any) -> bool:
    # name must match
    if not (isinstance(out_json, dict) and out_json.get("name")==dest_json.get("name")):
        return False
    # one protein + at least 1 ligand
    seqs = out_json.get("sequences")
    if not isinstance(seqs, list): return False
    prot = [el for el in seqs if isinstance(el, dict) and "protein" in el]
    lig  = [el for el in seqs if isinstance(el, dict) and "ligand"  in el]
    if len(prot)!=1 or len(lig)<1: return False
    protein = prot[0]["protein"]
    if not isinstance(protein, dict): return False
    # unpaired present; paired absent
    if "unpairedMsa" not in protein: return False
    if "pairedMsa" in protein: return False
    # no paired containers anywhere
    def has_paired(o):
        if isinstance(o, dict):
            for k,v in o.items():
                if k in DROP_EXACT_KEYS or k.lower().startswith("paired"):
                    return True
                if has_paired(v): return True
        elif isinstance(o, list):
            for v in o:
                if has_paired(v): return True
        return False
    if has_paired(out_json): return False
    # ligand equality check (structure must match)
    def first_lig(o):
        if isinstance(o, dict) and isinstance(o.get("sequences"), list):
            for el in o["sequences"]:
                if isinstance(el, dict) and "ligand" in el:
                    return el["ligand"]
        return None
    lig_dest = first_lig(dest_json) or find_ligand_value(dest_json)
    lig_out  = first_lig(out_json)  or find_ligand_value(out_json)
    if lig_dest is not None and lig_out is not None and lig_dest != lig_out:
        return False
    return True

def build_context_for_job(job_name: str, root: Path, donors_index: Dict[str, List[Path]]) -> Optional[JobContext]:
    dest_dir = root / job_name
    dest_input_json = dest_dir / "alphafold_input.json"
    if not dest_input_json.is_file():
        return None
    prefix = job_name.split("-")[0]
    donor_dir = find_donor(donors_index, prefix, exclude=job_name)
    if donor_dir is None:
        return None
    donor_with_msa_json = donor_dir / "output_msa" / "alphafold_input_with_msa.json"
    if not donor_with_msa_json.is_file():
        return None
    return JobContext(job_name, dest_dir, dest_input_json, donor_dir, donor_with_msa_json)

def make_single_chain_with_msa(ctx: JobContext) -> Path:
    dest = load_json(ctx.dest_input_json)
    donor = load_json(ctx.donor_with_msa_json)

    dest_seq = find_first_sequence(dest)
    keep_idx = choose_chain_index_by_sequence(donor, dest_seq)

    pruned = strip_paired_blocks(donor)
    pruned = prune_chain_dimension(pruned, keep_idx)

    # Replace ligand(s) from destination
    ligand_value = find_ligand_value(dest)
    out_json = override_ligand(pruned, ligand_value)

    # Set name to destination name
    if isinstance(dest, dict) and "name" in dest:
        out_json["name"] = dest["name"]

    # Drop optional keys
    for drop_key in ("bondedAtomPairs", "userCCD"):
        out_json.pop(drop_key, None)
    # Drop protein.modifications if present
    if isinstance(out_json.get("sequences"), list):
        for el in out_json["sequences"]:
            if isinstance(el, dict) and "protein" in el and isinstance(el["protein"], dict):
                el["protein"].pop("modifications", None)

    out_path = ctx.dest_dir / "output_msa" / "alphafold_input_with_msa.json"
    safe_write_json(out_path, out_json)
    return out_path

def finalize_or_fix_partial(out_path: Path, dest_json: Any) -> bool:
    # If a leftover .tmp exists and is valid, finalize it
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    if tmp.exists():
        try:
            tmp_json = load_json(tmp)
            if output_is_complete(tmp_json, dest_json):
                os.replace(tmp, out_path)
                return True
            else:
                tmp.unlink(missing_ok=True)
        except Exception:
            try: tmp.unlink(missing_ok=True)
            except Exception: pass
    # If final exists but invalid, let caller regenerate
    if out_path.exists():
        try:
            out_json = load_json(out_path)
            if output_is_complete(out_json, dest_json):
                return True
        except Exception:
            pass
    return False

def main():
    ap = argparse.ArgumentParser(description="Build single-chain alphafold_input_with_msa.json for all jobs, skipping completed ones.")
    ap.add_argument("--csv", type=Path, default=CSV_DEFAULT, help="CSV with one job name per line.")
    ap.add_argument("--root", type=Path, default=ROOT_DEFAULT, help="Root jobs/human_test_set directory.")
    ap.add_argument("--dry-run", action="store_true", help="Scan and report; do not write files.")
    ap.add_argument("--verbose", action="store_true", help="Print extra per-job details.")
    args = ap.parse_args()

    root = args.root
    if not root.exists():
        print(f"ERROR: root not found: {root}", file=sys.stderr)
        sys.exit(2)

    jobs = read_lines(args.csv)
    if not jobs:
        print(f"ERROR: no jobs found in {args.csv}", file=sys.stderr)
        sys.exit(2)

    donor_index = donors_by_prefix(root)
    summaries: List[str] = []

    for job in jobs:
        # skip header-like rows
        if job.lower() in {"input_folder_name", "job_name", "name"}:
            summaries.append(f"{job}: SKIP (header)")
            continue

        ctx = build_context_for_job(job, root, donor_index)
        if ctx is None:
            # Distinguish missing dest vs missing donor
            dest_input = root / job / "alphafold_input.json"
            if not dest_input.is_file():
                summaries.append(f"{job}: SKIP (no alphafold_input.json)")
            else:
                summaries.append(f"{job}: SKIP (no donor with MSA for prefix {job.split('-')[0]})")
            continue

        try:
            dest_json = load_json(ctx.dest_input_json)
        except Exception as e:
            summaries.append(f"{job}: ERROR (dest input unreadable: {e.__class__.__name__})")
            continue

        out_path = ctx.dest_dir / "output_msa" / "alphafold_input_with_msa.json"
        # If already complete, skip; else try to finalize from .tmp; else regenerate
        if out_path.exists():
            try:
                out_json = load_json(out_path)
                if output_is_complete(out_json, dest_json):
                    summaries.append(f"{job}: SKIP (already complete)")
                    continue
            except Exception:
                pass
            # Try finalize/cleanup partial
            if finalize_or_fix_partial(out_path, dest_json):
                summaries.append(f"{job}: FIXED (finalized partial)")
                continue
            # else regenerate below

        if args.dry_run:
            summaries.append(f"{job}: WOULD WRITE (donor={ctx.donor_dir.name}) -> {out_path}")
            continue

        try:
            out_path = make_single_chain_with_msa(ctx)
            # Validate immediately
            try:
                out_json = load_json(out_path)
                ok = output_is_complete(out_json, dest_json)
            except Exception:
                ok = False
            if ok:
                summaries.append(f"{job}: OK (donor={ctx.donor_dir.name})")
            else:
                summaries.append(f"{job}: ERROR (post-write validation failed)")
        except Exception as e:
            summaries.append(f"{job}: ERROR ({e.__class__.__name__}: {e})")

    # Print compact summary lines (no large content)
    print("\n".join(summaries))

if __name__ == "__main__":
    main()

