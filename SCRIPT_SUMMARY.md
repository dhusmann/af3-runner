# AlphaFold3 Pipeline Scripts Summary

## Complete Script List (Total: 28 scripts)

### Core Pipeline Scripts (7)
1. **batch_reuse_msa.py** - Identifies and reuses existing MSAs
2. **submit_msa_arrays.sh** - Distributes MSA jobs across partitions
3. **submit_msa_array.sh** - Runs individual MSA generation
4. **submit_dist.sh** - Manages GPU job distribution
5. **submit_gpu.sh** - Runs individual GPU inference
6. **launch_af3.sh** - Starts 48-hour automation
7. **af3_48hr_cycle.sh** - Handles periodic execution

### Monitoring Scripts (5)
8. **monitor_msa_arrays.sh** - Monitors MSA array jobs
9. **get_job_status.sh** - Reports job stages and seeds
10. **get_job_status_detailed.sh** - Advanced status with filters
11. **pipeline_status.sh** - Quick dashboard view
12. **pipeline_summary.sh** - Complete pipeline overview with sync status

### Output Sync Scripts (6)
13. **sync_organize_outputs.sh** - Syncs and organizes outputs locally
14. **rclone_to_gdrive.sh** - Submits SLURM job for Google Drive upload
15. **check_sync_status.sh** - Shows sync readiness
16. **sync_all.sh** - Complete sync workflow
17. **check_rclone_status.sh** - Monitor rclone SLURM jobs
18. **clean_output_dir.sh** - Cleans output directory

### Utility Scripts (4)
19. **cleanup_msa_tmp.sh** - Removes temporary MSA files
20. **pipeline_quickstart.sh** - Validates setup
21. **test_seed_detection.sh** - Tests seed completion detection
22. **test_rclone_quota.sh** - Tests quota error handling (optional)

### Helper Tools (Singles & maintenance) (6)
23. **tools/find_single_enzyme_ligand_jobs.py** - Identify single-protein–ligand jobs
24. **tools/clear_output_msa_jsons.py** - Remove `output_msa/*.json` for singles
25. **tools/restore_jobs_to_csv.py** - Merge verified singles into `folding_jobs.csv`
26. **tools/unixify_newlines.py** - Normalize CSVs to LF newlines
27. **tools/remove_outputs_from_list.py** - Delete `jobs/<name>/output` for a list
28. **tools/restore_jobs_from_list.py** - Add a list of jobs back to CSV

## Quick Setup

```bash
# Make all scripts executable
chmod +x *.sh *.py

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
