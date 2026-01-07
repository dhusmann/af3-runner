# Agent Notes (ChatGPT)

This repository includes helper scripts I added to streamline preparing single‑chain `alphafold_input_with_msa.json` files by borrowing the unpaired MSA from a related two‑chain donor job.

## What I Added

- `tools/msa_extract_chain1.py` — single‑job helper (development artifact) that:
  - Locates a donor job sharing the prefix before the first `-`.
  - Removes paired data and second chain, preserving only chain‑1 data.
  - Keeps `protein.unpairedMsa`, `protein.sequence`, and `protein.templates`.
  - Copies `name` and `ligand` from the destination `alphafold_input.json`.
  - Writes formatted JSON to `<dest>/output_msa/alphafold_input_with_msa.json`.

- `tools/msa_batch_extract_chain1.py` — batch script (recommended) that:
  - Processes all jobs listed in `folding_jobs_nsd2i.csv` under the given ROOT.
  - Skips already‑complete outputs and finalizes valid `*.tmp` files left by interrupted writes (atomic writes via `os.replace`).
  - Ensures only one protein, at least one ligand, unpaired MSA present, no paired content, destination `name` and ligand enforced.

- `tools/README.md` — documentation on purpose, assumptions, usage, validation, and troubleshooting.

## Key Assumptions / Behaviors

- Donor jobs contain two proteins (A/B) with per‑protein `unpairedMsa` and may include `pairedMsa` and complex‑level paired containers; all paired content is removed.
- Destination jobs are single‑chain; output must contain exactly one protein and at least one ligand.
- The ligand in outputs is always taken from the destination input (never copied from donor), and the top‑level `name` is set from the destination.
- Optional keys intentionally omitted from outputs: `protein.modifications`, top‑level `bondedAtomPairs`, and `userCCD`.

## Notable Functions (batch script)

- `strip_paired_blocks(obj)`: removes keys `pairedMsa`, `paired_msa(s)`, `complex_msa`, `interaction_msa` anywhere in the structure without touching `unpairedMsa`.
- `prune_chain_dimension(obj, keep_index)`: keeps only the first protein and all ligands in `sequences`; projects any per‑chain dicts (e.g., `seq: {A:..., B:...}`) to the kept chain.
- `output_is_complete(out, dest)`: validates one protein, ≥1 ligand, `unpairedMsa` present, no paired content, and `name`/ligand match destination.
- `safe_write_json(path, obj)`: atomic write with `*.tmp` then `os.replace`.

## Cleanup Guidance

The batch script writes atomically and should not leave partial files. If a run was interrupted, it may leave `alphafold_input_with_msa.json.tmp` under some `<job>/output_msa/` folders. Do not delete these manually unless confirmed; the batch script will auto‑finalize or remove invalid temps on the next run.

If you want me to perform any cleanup, please explicitly approve which files to remove.

---
Authored by ChatGPT (OpenAI).
