# Tools: Single‑Chain MSA Extraction

This folder contains helper scripts to construct per‑job `alphafold_input_with_msa.json` files for single‑chain runs by reusing the unpaired MSA from a related two‑chain donor job.

The scripts are designed to be safe on very large JSONs (avoid printing long alignment strings), resumeable after interruption, and idempotent (skip jobs that are already complete).

## Scripts

- `msa_batch_extract_chain1.py` (recommended)
  - Batch processes all jobs listed in a CSV (one job name per line).
  - For each destination job, finds a donor job sharing the prefix before the first `-` (e.g., `SMYD1-*`), extracts chain‑1 MSA, removes paired/second‑chain content, sets the correct `name` and `ligand`, and writes `output_msa/alphafold_input_with_msa.json` atomically.
  - Skips jobs that are already complete; finalizes any `*.tmp` left from an interrupted run.

- `msa_extract_chain1.py`
  - Single‑job helper used during development. It processes only the first valid job it finds in the CSV. Kept for reference.

## Fundamental Assumptions

- Destination (recipient) job:
  - Lives at: `<ROOT>/<job_name>/alphafold_input.json` (default ROOT is `/scratch/groups/ogozani/alphafold3/jobs/human_test_set`).
  - Is a single‑chain target: its output `alphafold_input_with_msa.json` must contain exactly one protein entry and at least one ligand entry in `sequences`.
  - Contains a `name` field (e.g., `SMYD1-SAH`) and a ligand object in `sequences` (or elsewhere) that should be used for the output.

- Donor job (provides MSA):
  - Resides under the same ROOT and shares the prefix before the first `-` with the destination (e.g., donor for `SMYD1-SAH` is some `SMYD1-*`).
  - Has `output_msa/alphafold_input_with_msa.json` produced from a two‑chain setup (or at least containing two protein entries and the per‑chain MSA fields).
  - Provides a per‑protein `unpairedMsa` field; `pairedMsa` may also be present (and will be removed).

- Chain mapping:
  - “Chain 1” means the first protein entry in donor data. The batch script heuristically confirms the match by comparing the destination’s first protein sequence against donor chains, and selects the better match. If uncertain, it keeps the first chain (index 0).

## What Gets Copied vs. Rewritten

- Copied/preserved from donor (for the selected chain only):
  - `protein.sequence`
  - `protein.unpairedMsa`
  - `protein.templates`
  - Any other single‑chain fields that are not paired‑specific

- Rewritten/derived from destination:
  - Top‑level `name` (set to destination’s `alphafold_input.json` name)
  - `sequences[].ligand` (copied from destination; never from donor)

- Removed/Omitted:
  - All paired content: the key `pairedMsa` and any containers with names like `paired_msas`, `paired_msa`, `complex_msa`, or `interaction_msa`.
  - The second chain (and any structures keyed by `{A,B}`, `{1,2}`, or `chain_1/chain_2` when applicable).
  - Optional keys: `modifications` (under protein), top‑level `bondedAtomPairs`, and `userCCD`.

## Idempotency, Resume, and Safety

- Atomic writes: files are written to a `*.tmp` first and moved into place on success.
- Resume: if a `*.tmp` exists and already validates, it is finalized automatically; otherwise it is discarded and rewritten.
- Skips completed jobs: the script validates that an existing output has exactly one protein, at least one ligand, includes `unpairedMsa`, contains no paired fields, and that the `name` and ligand match the destination.

## Usage

Default ROOT used by the scripts is `/scratch/groups/ogozani/alphafold3/jobs/human_test_set` and default CSV is `folding_jobs_nsd2i.csv` in the repo root.

- Dry‑run (no writes; prints one summary line per job):
```
python tools/msa_batch_extract_chain1.py \
  --csv folding_jobs_nsd2i.csv \
  --root /scratch/groups/ogozani/alphafold3/jobs/human_test_set \
  --dry-run
```

- Execute for all jobs (skips completed, finalizes partial):
```
python tools/msa_batch_extract_chain1.py \
  --csv folding_jobs_nsd2i.csv \
  --root /scratch/groups/ogozani/alphafold3/jobs/human_test_set
```

Notes:
- Python 3.8+ (std‑lib only). No third‑party deps.
- To force a rebuild for a single job, delete that job’s existing `output_msa/alphafold_input_with_msa.json` and re‑run.

## Validation (safe, no large prints)

Use `jq` to assert structure without dumping MSA strings:

- Name equals destination:
```
jq -r '.name' <dest/output_msa/alphafold_input_with_msa.json
jq -r '.name' <dest/alphafold_input.json
```

- Exactly one protein and at least one ligand:
```
jq '[.sequences[]|has("protein")]|map(select(.))|length' <dest/output_msa/alphafold_input_with_msa.json
jq '[.sequences[]|has("ligand")]|map(select(.))|length'  <dest/output_msa/alphafold_input_with_msa.json
```

- `unpairedMsa` present, `pairedMsa` absent (no values printed):
```
jq '..|objects|select(has("unpairedMsa"))|1' <dest/output_msa/alphafold_input_with_msa.json | head -n 1
jq '..|objects|select(has("pairedMsa"))|1'   <dest/output_msa/alphafold_input_with_msa.json | head -n 1
```

- No paired containers anywhere:
```
for k in complex_msa interaction_msa paired_msas paired_msa; do 
  grep -n -m 1 -i '"'$k'"' <dest/output_msa/alphafold_input_with_msa.json || echo "$k: absent"; 
done
```

## Troubleshooting

- “SKIP (no donor …)”: Add or generate a donor job with the same prefix (`PREFIX-*`) that already has `output_msa/alphafold_input_with_msa.json`.
- “ERROR (post‑write validation failed)”: Usually indicates the donor lacks `unpairedMsa` for the selected chain, or ligand/name mismatch. Inspect the donor’s protein keys (only the keys, not values):
```
jq -c '.sequences[]|select(has("protein"))|.protein|keys' <donor/output_msa/alphafold_input_with_msa.json | head -n 3
```
- Chain mismatch: If the heuristic picked the wrong donor chain, ensure the destination’s protein sequence in `alphafold_input.json` matches the intended donor chain; otherwise, point the destination job to a closer donor.

---
If you want different donor selection logic or a `--force` flag, let me know and I’ll add it.

