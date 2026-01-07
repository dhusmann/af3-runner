# AlphaFold3 Pipeline Scripts Summary

This repo contains scripts for:
- Creating AlphaFold3 job folders + `alphafold_input.json`
- Running a two-stage pipeline on Sherlock/Slurm (MSA → GPU inference), with MSA reuse
- Monitoring job progress
- Organizing outputs + archiving MSAs + compressing seeds + syncing to Google Drive

## Key Entrypoints

- **Create jobs**: `./createAF3query.sh --help`
- **Run pipeline (automated cycles)**: `./launch_af3.sh`
- **Sync + backup (compute-node, parallel)**: `./sync_all.sh`

## Script List (Total: 45)

### Job Creation (4)
- `createAF3query.sh` - Unified job creation (PTMs incl. ALL/EACH; CCD + SMILES ligands)
- `createHTS-AF3query.sh` - Wrapper: defaults to `jobs/human_test_set/` + `folding_jobs.csv`
- `createAF3query_withSMILES.sh` - Wrapper: defaults to `jobs/human_test_set/` + `folding_jobs_nsd2i.csv`
- `createBatchAF3queries.sh` - Batch helper that repeatedly calls `createAF3query.sh`

### Core Pipeline (7)
- `batch_reuse_msa.py` - MSA reuse + triage (`msa_array_jobs.csv`, `waiting_for_msa.csv`)
- `submit_msa_arrays.sh` - Split + submit MSA arrays across partitions
- `submit_msa_array.sh` - MSA-only worker (`--norun_inference`)
- `submit_dist.sh` - GPU submitter (runs `batch_reuse_msa.py`, drains waiting list)
- `submit_gpu.sh` - GPU worker (multi-chain inference-only; singles run full pipeline)
- `launch_af3.sh` - Starts the automated cycling pipeline
- `af3_48hr_cycle.sh` - SLURM cycle runner that re-submits itself

### Monitoring (5)
- `monitor_msa_arrays.sh` - MSA array monitor (`squeue`/`sacct` + quick hints)
- `get_job_status.sh` - Stage scan + seed completion heuristic
- `get_job_status_detailed.sh` - Adds filters + CSV export
- `pipeline_status.sh` - Compact dashboard
- `pipeline_summary.sh` - Expanded dashboard (includes sync readiness)

### Output Organization + Sync (17)
- `sync_all.sh` - One-command orchestrator (local processing + optional gdrive sync)
- `sync_organize_outputs.sh` - Submits compute-node jobs (rsync, MSA archive, seed pack)
- `sync_parallel.conf` - Parallelism + resource config for sync
- `sync_job_discovery.sh` - Discovers output subdirs and chunks work for arrays
- `sync_organize_rsync.sbatch` - SLURM array wrapper for rsync chunks
- `sync_organize_rsync.sh` - Parallel rsync worker + destination routing
- `archive_msa_data.sbatch` - SLURM array wrapper (one job-group per task)
- `archive_msa_data.sh` - Deduplicating MSA archiver → `output/msa/*.tar.gz`
- `compress_seeds_array.sbatch` - SLURM array wrapper (one output dir per task)
- `compress_seeds.sh` - Creates `seeds.tar.gz` per run dir (keeps originals)
- `compress_seeds.sbatch` - Older “one-job does everything” seed packing
- `pack_seeds_human_test_set.sbatch` - Human-test-set-only seed packing
- `check_sync_status.sh` - Reports what’s ready to sync (sizes, counts)
- `clean_output_dir.sh` - Interactive reset of `output/`
- `rclone_to_gdrive.sh` - Submits gdrive sync job (also backs up scripts/tools/analysis)
- `rclone_retry.sh` - Retry wrapper for gdrive sync (quota handling)
- `check_rclone_status.sh` - Monitor gdrive sync jobs + recent logs
- `rclone_msa_to_gdrive.sh` - Legacy MSA-only rclone job (generally superseded)

### Utilities (2)
- `cleanup_msa_tmp.sh` - Removes `msa_array_jobs_part*.tmp` after arrays finish
- `pipeline_quickstart.sh` - Validates setup (paths + required scripts)

### Tools (10)
- `tools/find_single_enzyme_ligand_jobs.py` - Verify singles list (1 protein + 1 ligand)
- `tools/clear_output_msa_jsons.py` - Remove `output_msa/*.json` for verified singles
- `tools/restore_jobs_to_csv.py` - Merge verified singles into `folding_jobs.csv`
- `tools/unixify_newlines.py` - Convert pipeline CSVs to LF newlines
- `tools/remove_outputs_from_list.py` - Remove `output/` dirs for a job list
- `tools/restore_jobs_from_list.py` - Add a list back to `folding_jobs.csv`
- `tools/msa_batch_extract_chain1.py` - Build single-chain `alphafold_input_with_msa.json` from donor jobs
- `tools/msa_extract_chain1.py` - Single-job helper (development artifact)
- `tools/move_msa_data.sh` - One-off helper to move `*_data.json` into `output_msa/`
- `tools/remove_empty_outputs.sh` - One-off helper to remove empty `output/` dirs

## Quick Setup

```bash
# Make all scripts executable
chmod +x *.sh *.sbatch *.py tools/*.sh tools/*.py

# Verify setup
./pipeline_quickstart.sh
```

## Typical Workflow

1. **Start pipeline**: `./launch_af3.sh`
2. **Monitor progress**: `./pipeline_status.sh`
3. **Check job details**: `./get_job_status.sh`
4. **Sync outputs**: `./sync_all.sh`

## Script Categories by Function

### Starting Jobs
- `launch_af3.sh` → `af3_48hr_cycle.sh` → `submit_dist.sh`
- `submit_msa_arrays.sh` → `submit_msa_array.sh`
- `submit_gpu.sh`

### Creating Jobs
- `createAF3query.sh` (unified)
- `createHTS-AF3query.sh` (HTS defaults wrapper)
- `createAF3query_withSMILES.sh` (SMILES defaults wrapper)

## Updates (2025-09-17)
- `submit_dist.sh` (updated):
  - Submits single-protein–ligand jobs without requiring `output_msa/`.
  - Strips CR from CSV lines to avoid silent directory mismatches.
  - Applies higher CPU/time for singles.
- `submit_gpu.sh` (updated):
  - Runs AF3 data pipeline for singles from `alphafold_input.json` (no `--norun_data_pipeline`).
  - Keeps GPU-only for multi-chain and prefers `output_msa/alphafold_input_with_msa.json`.
- `batch_reuse_msa.py` (updated):
  - Excludes single-protein–ligand jobs from MSA reuse/copy; writes CSVs with LF newlines.

### Monitoring Progress
- `pipeline_status.sh` - Overall view
- `monitor_msa_arrays.sh` - MSA jobs
- `get_job_status.sh` - Individual jobs
- `get_job_status_detailed.sh` - Advanced queries

### Managing Outputs
- `check_sync_status.sh` - Pre-sync check
- `sync_organize_outputs.sh` - Local organization
- `rclone_to_gdrive.sh` - Cloud upload
- `sync_all.sh` - Complete workflow

### Maintenance
- `cleanup_msa_tmp.sh` - Clean MSA temps
- `clean_output_dir.sh` - Reset output directory
- `pipeline_quickstart.sh` - Verify installation
