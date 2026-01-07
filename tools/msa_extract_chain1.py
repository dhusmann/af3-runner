#!/usr/bin/env python3
from __future__ import annotations
"""
Create alphafold_input_with_msa.json for single-chain jobs by borrowing
chain-1 MSA from a related two-chain job sharing the same prefix.

Behavior:
- Reads the first job_name from a CSV (one per line).
- Destination: jobs/human_test_set/{job_name}/alphafold_input.json
- Donor: another directory under jobs/human_test_set/ that starts with the
  {prefix} before the first hyphen (e.g., SMYD1-...), and contains
  output_msa/alphafold_input_with_msa.json (two-chain source).
- Produces: jobs/human_test_set/{job_name}/output_msa/alphafold_input_with_msa.json

Rules:
- Only include chain 1 (first chain) from the donor; remove any second chain and any paired MSA.
- The ligand in the output must be taken from the destination alphafold_input.json,
  not from the donor with_msa JSON.
- Write pretty-printed JSON with indent=2.

Notes:
- The code uses schema introspection and heuristics to find chain containers and paired MSA blocks
  without assuming exact key names. It aims to be conservative: if it isn't confident something
  is chain-like, it leaves it as-is.
"""

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


ROOT_DEFAULT = Path("jobs/human_test_set")


def read_lines(path: Path) -> List[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def load_json(path: Path) -> Any:
    with path.open("r") as f:
        return json.load(f)


def dump_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
        f.write("\n")


def is_aa_seq(s: str) -> bool:
    return bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", s)) and len(s) >= 30


def find_first_sequence(obj: Any) -> Optional[str]:
    best: Optional[str] = None

    def visit(x: Any):
        nonlocal best
        if isinstance(x, dict):
            for v in x.values():
                visit(v)
        elif isinstance(x, list):
            for v in x:
                visit(v)
        elif isinstance(x, str):
            if is_aa_seq(x):
                if best is None or len(x) > len(best):
                    best = x

    visit(obj)
    return best


def find_ligand_value(obj: Any) -> Optional[Any]:
    q: List[Any] = [obj]
    while q:
        cur = q.pop(0)
        if isinstance(cur, dict):
            if "ligand" in cur:
                return cur.get("ligand")
            if "ligands" in cur:
                return cur.get("ligands")
            q.extend(cur.values())
        elif isinstance(cur, list):
            q.extend(cur)
    return None


def listdirs(p: Path) -> List[Path]:
    return [x for x in p.iterdir() if x.is_dir()]


def find_donor_with_prefix(root: Path, prefix: str, exclude: str) -> Optional[Path]:
    candidates = []
    for d in listdirs(root):
        name = d.name
        if name == exclude:
            continue
        if not name.startswith(prefix + "-"):
            continue
        donor_json = d / "output_msa" / "alphafold_input_with_msa.json"
        if donor_json.is_file():
            candidates.append(d)
    candidates.sort(key=lambda p: p.name)
    return candidates[0] if candidates else None


def looks_chain_like_element(x: Any) -> bool:
    if not isinstance(x, dict):
        return False
    keys = set(x.keys())
    hints = {"msa", "templates", "sequence", "seq", "chain_id", "alignments"}
    return bool(keys & hints)


def is_chain_like_list(lst: List[Any]) -> bool:
    if len(lst) < 2:
        return False
    dicts = [x for x in lst if isinstance(x, dict)]
    if len(dicts) < max(2, len(lst) // 2):
        return False
    return any(looks_chain_like_element(x) for x in dicts)


PAIR_CONTAINER_KEYS = {"paired_msa", "paired_msas", "complex_msa", "interaction_msa"}
DROP_EXACT_KEYS = {"pairedMsa"} | PAIR_CONTAINER_KEYS


def strip_paired_blocks(d: Any) -> Any:
    if isinstance(d, dict):
        out = {}
        for k, v in d.items():
            # Drop only explicit paired keys/containers; DO NOT drop 'unpairedMsa'
            if k in DROP_EXACT_KEYS or k.lower().startswith("paired"):
                continue
            out[k] = strip_paired_blocks(v)
        return out
    elif isinstance(d, list):
        return [strip_paired_blocks(x) for x in d]
    return d


def choose_chain_index_by_sequence(donor: Any, dest_seq: Optional[str]) -> int:
    if not dest_seq:
        return 0
    candidates: List[List[Any]] = []
    if isinstance(donor, dict):
        for k in ("chains", "sequences", "entities", "monomers", "inputs", "targets"):
            v = donor.get(k)
            if isinstance(v, list) and is_chain_like_list(v):
                candidates.append(v)
    if not candidates:
        def scan(obj: Any):
            if isinstance(obj, dict):
                for vv in obj.values():
                    scan(vv)
            elif isinstance(obj, list):
                if is_chain_like_list(obj):
                    candidates.append(obj)
                else:
                    for vv in obj:
                        scan(vv)
        scan(donor)

    def seq_from_elem(elem: Any) -> Optional[str]:
        if isinstance(elem, dict):
            for key in ("sequence", "seq", "query_sequence", "target_sequence"):
                s = elem.get(key)
                if isinstance(s, str) and is_aa_seq(s):
                    return s
        return None

    best_idx = 0
    best_score = -1
    for lst in candidates:
        for idx in range(min(2, len(lst))):
            s = seq_from_elem(lst[idx])
            score = 0
            if s:
                if dest_seq in s or s in dest_seq:
                    score = 3
                else:
                    n = min(len(s), len(dest_seq))
                    if n:
                        matches = sum(1 for a, b in zip(s[:n], dest_seq[:n]) if a == b)
                        score = int(100 * matches / max(1, n))
            if score > best_score:
                best_score, best_idx = score, idx
    return best_idx


def prune_chain_dimension(obj: Any, keep_index: int, chain_keys: Tuple[str, str] = ("A", "B")) -> Any:
    if isinstance(obj, dict):
        keys = set(obj.keys())
        ab_like = keys.issuperset(set(chain_keys))
        one_two_like = keys.issuperset({"1", "2"}) or keys.issuperset({"chain_1", "chain_2"})
        if ab_like or one_two_like:
            key_keep = chain_keys[keep_index] if ab_like else ("1" if keep_index == 0 else "2")
            if one_two_like and ("chain_1" in keys or "chain_2" in keys):
                key_keep = "chain_1" if keep_index == 0 else "chain_2"
            if key_keep in obj:
                return prune_chain_dimension(obj[key_keep], keep_index, chain_keys)

        out: Dict[str, Any] = {}
        for k, v in obj.items():
            if k in DROP_EXACT_KEYS or k.lower().startswith("paired"):
                continue
            if k == "seq" and isinstance(v, dict):
                for cand in (chain_keys[keep_index], str(keep_index + 1)):
                    if cand in v and isinstance(v[cand], (str, dict, list)):
                        out[k] = prune_chain_dimension(v[cand], keep_index, chain_keys)
                        break
                else:
                    out[k] = prune_chain_dimension(v, keep_index, chain_keys)
                continue
            out[k] = prune_chain_dimension(v, keep_index, chain_keys)
        return out

    if isinstance(obj, list):
        # Special-case: sequences arrays of {protein:{...}} and {ligand:{...}}
        if obj and all(isinstance(x, dict) for x in obj) and any(("protein" in x or "ligand" in x) for x in obj):
            out_list: List[Any] = []
            prot_seen = 0
            for el in obj:
                if "ligand" in el:
                    out_list.append(prune_chain_dimension(el, keep_index, chain_keys))
                elif "protein" in el:
                    if prot_seen == keep_index:
                        out_list.append(prune_chain_dimension(el, keep_index, chain_keys))
                    prot_seen += 1
            return out_list
        if is_chain_like_list(obj) and len(obj) > keep_index:
            return prune_chain_dimension(obj[keep_index], keep_index, chain_keys)
        return [prune_chain_dimension(x, keep_index, chain_keys) for x in obj]

    return obj


def override_ligand(out_json: Any, ligand_value: Any) -> Any:
    if ligand_value is None:
        return out_json
    def walk(x: Any) -> Any:
        if isinstance(x, dict):
            y: Dict[str, Any] = {}
            for k, v in x.items():
                if k == "ligand":
                    y[k] = ligand_value
                elif k == "ligands":
                    y[k] = ligand_value if isinstance(ligand_value, list) else [ligand_value]
                else:
                    y[k] = walk(v)
            return y
        if isinstance(x, list):
            return [walk(v) for v in x]
        return x
    return walk(out_json)


@dataclass
class JobContext:
    job_name: str
    dest_dir: Path
    dest_input_json: Path
    donor_dir: Path
    donor_with_msa_json: Path


def build_context_for_first_job(csv_path: Path, root: Path) -> JobContext:
    jobs = read_lines(csv_path)
    if not jobs:
        raise SystemExit(f"No jobs found in {csv_path}")
    # Choose the first job that actually has a destination input JSON present
    job_name = None
    for j in jobs:
        cand = root / j / "alphafold_input.json"
        if cand.is_file():
            job_name = j
            dest_input_json = cand
            dest_dir = dest_input_json.parent
            break
    if job_name is None:
        raise SystemExit(f"No valid destination alphafold_input.json found under {root} for any job in {csv_path}")
    prefix = job_name.split("-")[0]
    donor_dir = find_donor_with_prefix(root, prefix, exclude=job_name)
    if donor_dir is None:
        raise SystemExit(f"No donor found under {root} with prefix {prefix}-*")
    donor_with_msa_json = donor_dir / "output_msa" / "alphafold_input_with_msa.json"
    return JobContext(job_name, dest_dir, dest_input_json, donor_dir, donor_with_msa_json)


def make_single_chain_with_msa(ctx: JobContext, dry_run: bool = False) -> Path:
    dest = load_json(ctx.dest_input_json)
    donor = load_json(ctx.donor_with_msa_json)

    dest_seq = find_first_sequence(dest)
    keep_idx = choose_chain_index_by_sequence(donor, dest_seq)

    pruned = strip_paired_blocks(donor)
    pruned = prune_chain_dimension(pruned, keep_idx)

    ligand_value = find_ligand_value(dest)
    out_json = override_ligand(pruned, ligand_value)
    # Ensure output name matches destination job's name
    if isinstance(dest, dict) and "name" in dest:
        out_json["name"] = dest["name"]

    # Drop optional keys per requirement
    for drop_key in ("bondedAtomPairs", "userCCD"):
        if drop_key in out_json:
            out_json.pop(drop_key, None)
    # Drop protein.modifications if present
    try:
        if isinstance(out_json.get("sequences"), list):
            for el in out_json["sequences"]:
                if isinstance(el, dict) and "protein" in el and isinstance(el["protein"], dict):
                    el["protein"].pop("modifications", None)
    except Exception:
        pass

    out_path = ctx.dest_dir / "output_msa" / "alphafold_input_with_msa.json"

    if dry_run:
        print("[DRY-RUN] Would write:", out_path)
    else:
        dump_json(out_path, out_json)
        print("Wrote:", out_path)

    return out_path


def main():
    ap = argparse.ArgumentParser(description="Create single-chain alphafold_input_with_msa.json from donor two-chain jobs")
    ap.add_argument("--csv", default="folding_jobs_nsd2i.csv", type=Path, help="CSV with one job_name per line (default: folding_jobs_nsd2i.csv)")
    ap.add_argument("--root", default=str(ROOT_DEFAULT), type=Path, help="Root folder containing jobs/human_test_set (default: jobs/human_test_set)")
    ap.add_argument("--dry-run", action="store_true", help="Do not write output, just report actions")
    ap.add_argument("--all", action="store_true", help="Process all jobs in the CSV instead of just the first valid one")
    args = ap.parse_args()

    root = args.root
    if root.name != "human_test_set":
        if root.name == "jobs" and (root / "human_test_set").is_dir():
            root = root / "human_test_set"

    if not args.all:
        ctx = build_context_for_first_job(args.csv, root)
        print(f"Destination job: {ctx.job_name}")
        print(f"  dest input:   {ctx.dest_input_json}")
        print(f"Donor job:      {ctx.donor_dir.name}")
        print(f"  donor msa:    {ctx.donor_with_msa_json}")
        make_single_chain_with_msa(ctx, dry_run=args.dry_run)
        return

    # Batch mode: process all jobs listed
    jobs = read_lines(args.csv)
    summaries: List[str] = []
    for job_name in jobs:
        dest_input_json = root / job_name / "alphafold_input.json"
        if not dest_input_json.is_file():
            summaries.append(f"{job_name}: SKIP (no alphafold_input.json)")
            continue
        prefix = job_name.split("-")[0]
        donor_dir = find_donor_with_prefix(root, prefix, exclude=job_name)
        if donor_dir is None:
            summaries.append(f"{job_name}: SKIP (no donor)")
            continue
        ctx = JobContext(job_name, dest_input_json.parent, dest_input_json, donor_dir, donor_dir / "output_msa" / "alphafold_input_with_msa.json")
        try:
            out_path = make_single_chain_with_msa(ctx, dry_run=args.dry_run)
            summaries.append(f"{job_name}: OK (donor={donor_dir.name}) -> {out_path}")
        except Exception as e:
            summaries.append(f"{job_name}: ERROR ({e.__class__.__name__}: {e})")

    print("\n".join(summaries))


if __name__ == "__main__":
    main()
