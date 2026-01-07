# AlphaFold3 High-Throughput Pipeline

A comprehensive, automated pipeline for running AlphaFold3 predictions at scale on HPC clusters with intelligent MSA reuse and multi-partition job distribution.

## Changelog

- 2025-09-17
  - Added first-class support for single-protein–ligand jobs (ENZYME-SAM/SAH): these now run AF3’s data pipeline during the GPU job instead of reusing MSAs.
  - Updated `submit_dist.sh` to detect singles, submit them without requiring `output_msa/`, strip CRLF newlines when reading CSVs, and request more CPU/time for singles.
  - Updated `submit_gpu.sh` to run AF3 with the full data pipeline for singles from `alphafold_input.json` (keeps GPU-only for multi-chain).
  - Updated `batch_reuse_msa.py` to skip MSA reuse/copy for singles and write CSVs with Unix LF.
  - Added helper tools in `tools/` for discovering singles, cleaning existing MSA/outputs, restoring jobs to CSV, and normalizing CSV newlines.
  - Documented an operational caveat: non-empty `output/` may falsely mark a job as complete; provided cleanup tools and noted a future refinement.

## Overview

This pipeline processes large batches of AlphaFold3 predictions by:
1. **Reusing existing MSAs** from similar protein pairs to avoid redundant computation
2. **Distributing MSA generation** across multiple SLURM partitions using array jobs
3. **Managing GPU inference jobs** within cluster quotas
4. **Automating the workflow** with periodic job submission over 48 hours

## Pipeline Architecture

```
folding_jobs.csv
    ↓
[CYCLE 1 ONLY]:
batch_reuse_msa.py ─────┬─→ Copies existing MSAs
    ↓                    ├─→ msa_array_jobs.csv (new MSAs needed)
    ↓                    └─→ waiting_for_msa.csv (dependent jobs)
    ↓
submit_msa_arrays.sh ───────→ Array jobs on normal + hns partitions
    ↓
[ALL CYCLES]:
submit_dist.sh ─────────────→ GPU inference (up to 48 concurrent)
    ↑                           (includes batch_reuse_msa.py)
    └── Every 2 hours via af3_48hr_cycle.sh (24 cycles)
```

## Components

### Core Scripts

1. **`batch_reuse_msa.py`** - MSA reuse and job triage
   - Scans all jobs to identify which can reuse existing MSAs
   - Creates `msa_array_jobs.csv` for jobs needing fresh MSAs
   - Creates `waiting_for_msa.csv` for dependent jobs

2. **`submit_msa_arrays.sh`** - MSA array job distribution
   - Splits jobs between `normal` and `hns` partitions
   - Handles jobs exceeding SLURM array limits (1000)
   - Creates persistent `.tmp` files for each partition

3. **`submit_msa_array.sh`** - Individual MSA job execution
   - Runs MSA generation for a single job
   - Handles jobs in multiple directory locations
   - Updates central log file

4. **`submit_dist.sh`** - GPU job distribution
   - Integrates MSA reuse by running `batch_reuse_msa.py` first
   - Processes both main queue and waiting jobs
   - Respects GPU quota limits (default: 48)
   - Removes completed jobs from CSVs

5. **`submit_gpu.sh`** - Individual GPU inference job
   - Runs AlphaFold3 inference on GPU
   - Uses augmented JSON from MSA stage
   - Logs completion with timing information

### Automation Scripts

6. **`launch_af3.sh`** - Pipeline launcher
   - Initiates the 48-hour cycling process
   - Submits first cycle job

7. **`af3_48hr_cycle.sh`** - Periodic execution
   - Runs `submit_dist.sh` every 2 hours
   - Self-submits next cycle using SLURM dependencies
   - Sends completion email after 24 cycles

### Monitoring Scripts

8. **`monitor_msa_arrays.sh`** - MSA job monitoring
   - Shows job counts by partition
   - Displays recent failures
   - Provides useful SLURM commands

9. **`get_job_status.sh`** - Comprehensive job status report
   - Checks all jobs or specific jobs from CSV
   - Reports job stages (1: need MSA, 2: need GPU, 3: complete)
   - Lists completed seeds for Stage 3 jobs
   - Provides summary statistics

10. **`get_job_status_detailed.sh`** - Extended status reporting
    - All features of get_job_status.sh plus:
    - Export results to CSV
    - Filter by stage or seed
    - Show completion timestamps
    - Display job directory paths

11. **`pipeline_status.sh`** - Quick pipeline dashboard
    - Shows if automation is running and current cycle
    - Displays MSA and GPU job counts
    - Provides quick progress estimate
    - Shows recent completions

### Utility Scripts

19. **`cleanup_msa_tmp.sh`** - Temporary file cleanup
    - Removes partition `.tmp` files
    - Warns if jobs are still running

20. **`pipeline_quickstart.sh`** - Quick setup and validation
    - Verifies all scripts are present
    - Checks directory structure
    - Makes scripts executable

21. **`test_seed_detection.sh`** - Debug seed detection
    - Tests seed completion checking
    - Helps troubleshoot Stage 3 jobs

22. **`test_rclone_quota.sh`** - Test quota handling
    - Verifies quota error detection
    - Optional testing script

### Output Sync Scripts

13. **`sync_organize_outputs.sh`** - Sync and organize outputs
    - Copies outputs from `jobs/` to `output/`
    - Excludes large `*_data.json` files
    - Organizes by project (arabidopsis, ecoli, yeast, etc.)
    - Options for quiet/verbose output

14. **`rclone_to_gdrive.sh`** - Submit Google Drive sync job
    - Submits SLURM job to hns partition
    - Handles Google Drive quota limits
    - Auto-reschedules if quota exceeded
    - Sends email notifications

15. **`check_sync_status.sh`** - Check sync readiness
    - Shows which jobs have outputs ready
    - Calculates total size to sync
    - Reports excluded file sizes

16. **`sync_all.sh`** - Complete sync workflow
    - Combines all sync steps
    - Includes safety prompts
    - Handles both local and cloud sync

17. **`check_rclone_status.sh`** - Monitor rclone jobs
    - Shows active SLURM sync jobs
    - Reports completed/failed syncs
    - Checks for quota errors

18. **`clean_output_dir.sh`** - Clean output directory
    - Shows current usage
    - Safely removes output directory


## Updates (2025-09-17) — Single-Protein–Ligand Support

Overview
- Added first-class handling for single-protein–ligand jobs named `ENZYME-SAM` or `ENZYME-SAH` (exactly one protein and one ligand in `alphafold_input.json`).
- These jobs now run the AF3 data pipeline (MSA building) inside the GPU job; they do not reuse/copy MSAs from other jobs.
- Multi-chain enzyme–substrate–ligand jobs are unchanged and continue to use the two-stage CPU (MSA) → GPU flow with reuse.

Scripts changed
- `submit_dist.sh`
  - Detects single-protein–ligand jobs and submits them even when `output_msa/` is empty.
  - Normalizes Windows newlines (CRLF) in CSV lines when reading to avoid silent skips.
  - Singles request more resources when submitted (see Resources below).
- `submit_gpu.sh`
  - Singles: runs AF3 with the full data pipeline (no `--norun_data_pipeline`) using the bound job JSON (`alphafold_input.json`).
  - Multi-chain: unchanged; still prefers `output_msa/alphafold_input_with_msa.json` and uses `--norun_data_pipeline`.
  - Binds the job directory into the container for singles and sets `OMP_NUM_THREADS` from `SLURM_CPUS_PER_TASK`.
- `batch_reuse_msa.py`
  - Skips MSA reuse for single-protein–ligand jobs so they always build fresh MSAs inside AF3.
  - Writes CSVs (`msa_array_jobs.csv`, `waiting_for_msa.csv`) with Unix LF newlines.

New helper tools (under `tools/`)
- `find_single_enzyme_ligand_jobs.py`: Produce `analysis/single_enzyme_ligand_jobs.list/.tsv` by verifying 1 protein + 1 ligand.
- `clear_output_msa_jsons.py`: Remove `output_msa/*.json` for the verified singles.
- `restore_jobs_to_csv.py`: Add verified singles back into `folding_jobs.csv` (deduped, LF newlines).
- `unixify_newlines.py`: Convert `folding_jobs.csv`, `waiting_for_msa.csv`, `msa_array_jobs.csv` to LF to avoid CRLF issues.
- `remove_outputs_from_list.py`: Remove `jobs/human_test_set/<name>/output` for each name in a list (e.g., snapshots) to force resubmission when partial outputs exist.
- `restore_jobs_from_list.py`: Merge arbitrary job lists back into `folding_jobs.csv`.

Single-protein–ligand workflow
- Stage 1 (MSA arrays): Skipped for singles (no reuse/copy).
- Stage 2 (GPU): Submitted even if `output_msa/` is empty; AF3 runs its data pipeline to build MSAs on allocated CPUs.
- Stage 3 (complete): Same outputs as before under `jobs/<name>/output`.

Resources for singles
- Singles run data pipeline + inference in the GPU job, so they request:
  - `cpus-per-task`: 12
  - `time`: 08:00:00
- Multi-chain jobs keep previous defaults.

Operational caveat and workaround
- The current “GPU complete” check in `submit_dist.sh` considers a job complete if `jobs/<name>/output` exists and is non-empty. Partial leftovers can cause a job to be skipped.
- Workaround: use `tools/remove_outputs_from_list.py` with a job list snapshot to remove `output/` for affected singles, then `tools/restore_jobs_from_list.py` to add them back to `folding_jobs.csv`. A future refinement may switch the completeness test to look for AF3 completion markers (e.g., `summary_confidences.json`).

## Directory Structure

```
/scratch/groups/ogozani/alphafold3/
├── folding_jobs.csv              # Main job list
├── msa_array_jobs.csv            # Jobs needing MSA generation
├── waiting_for_msa.csv           # Jobs waiting for MSAs
├── logged_folding_jobs.csv       # Comprehensive job log
├── msa_array_jobs_part1.tmp     # Normal partition jobs
├── msa_array_jobs_part2.tmp     # HNS partition jobs
├── jobs/
│   ├── {job_name}/
│   │   ├── alphafold_input.json
│   │   ├── output_msa/          # MSA results
│   │   │   └── *_data.json
│   │   └── output/              # GPU inference results
│   └── human_test_set/
│       └── {job_name}/          # Alternative location
└── logs/
    ├── {job_id}_{array_index}_MSA.out/err
    └── {job_id}_GPU-only.out/err
```

## Input Requirements

### folding_jobs.csv
```csv
input_folder_name
SETD6-CALM1_noM-SAH
SETD6-CALM1_noM_K115me1-SAH
EEF1AKMT1-ACTB-SAH
```

### Job Naming Convention
- Format: `{entry1}-{entry2}[_PTM][-{cofactor}]`
- PTM markers: `_K`, `_H`, `_Q` followed by position and modification
- Example: `SETD6-CALM1_noM_K115me1-SAH`

## Usage

### One-Time Setup
```bash
# Make all scripts executable
chmod +x *.sh *.py

# Verify python module is available
ml python/3.9.0

# Optional: Test the monitoring scripts
./pipeline_status.sh              # Quick overview
./get_job_status.sh -h           # See all options
```

### Running the Complete Pipeline

#### Option 1: Automated 48-Hour Run (Recommended)
```bash
# Start the automated pipeline (runs for 48 hours)
./launch_af3.sh
```

**What happens automatically:**
1. **Cycle 1** (immediately):
   - Runs `batch_reuse_msa.py` to analyze all jobs
   - Submits MSA array jobs if needed via `submit_msa_arrays.sh`
   - May skip GPU submission if many MSA jobs are queued
2. **Cycles 2-24** (every 2 hours):
   - Runs `batch_reuse_msa.py` to check for newly available MSAs
   - Processes jobs from `waiting_for_msa.csv`
   - Submits GPU jobs up to quota limit
3. **After 48 hours**:
   - Sends completion email
   - All logs available for review

```bash
# Monitor progress
squeue -u $USER
tail -f af3_cycle_*.out
```

#### Option 2: Manual Step-by-Step
```bash
# 1. Analyze and reuse existing MSAs
python batch_reuse_msa.py

# 2. Submit MSA generation array jobs
./submit_msa_arrays.sh msa_array_jobs.csv

# 3. Monitor MSA progress
./monitor_msa_arrays.sh

# 4. Submit GPU jobs (run periodically)
./submit_dist.sh

# 5. Clean up temporary files when done
./cleanup_msa_tmp.sh
```

### Monitoring

```bash
# Quick pipeline overview
./pipeline_status.sh

# Detailed job status
./get_job_status.sh                           # All jobs
./get_job_status.sh -f specific_jobs.csv      # Specific jobs
./get_job_status.sh -s                        # Summary only
./get_job_status.sh -v                        # Verbose with seed details

# Advanced status reporting
./get_job_status_detailed.sh -stage 2         # Show only Stage 2 jobs
./get_job_status_detailed.sh -e report.csv    # Export to CSV
./get_job_status_detailed.sh -seed 0          # Jobs with seed 0 complete
./get_job_status_detailed.sh -p               # Show full paths

# Monitor MSA arrays specifically
./monitor_msa_arrays.sh

# Check GPU jobs
squeue -u $USER -p gpu

# View logs
tail -f logs/*_MSA.out
tail -f *_GPU-only.out

# Check pipeline progress
tail -f logged_folding_jobs.csv
```

## Job Stages

The pipeline tracks each job through three stages:

1. **Stage 1 - Need MSA**: Job has `alphafold_input.json` but no MSA output yet
2. **Stage 2 - Need GPU**: MSA generation complete, ready for GPU inference
3. **Stage 3 - Complete**: GPU inference finished, structure predictions available

Use `get_job_status.sh` to see which stage each job is in.

## Key Features

### Intelligent MSA Reuse
- Identifies protein pairs with existing MSAs
- Copies and merges MSA data to avoid recomputation
- Handles post-translational modifications intelligently

### Multi-Partition Distribution
- Splits MSA jobs between `normal` and `hns` partitions
- Utilizes separate quotas for maximum throughput
- Handles arrays exceeding 1000-job limit

### Automated Workflow
- Processes `waiting_for_msa.csv` jobs as MSAs complete
- Continuously submits GPU jobs within quota limits
- **Cycle 1**: Runs initial MSA analysis and submits array jobs
- **Cycles 2-24**: Processes GPU jobs every 2 hours
- Runs for 48 hours total with email notification

### Robust Error Handling
- Validates job directories and input files
- Handles missing directories gracefully
- Logs all operations for debugging

## Configuration

### Modifying Resource Allocations

**MSA Jobs** (`submit_msa_array.sh`):
```bash
#SBATCH --ntasks-per-node=32
#SBATCH --mem=128GB
#SBATCH --time=8:00:00
```

**GPU Jobs** (`submit_gpu.sh`):
```bash
#SBATCH --mem=64GB
#SBATCH --gpus=1
#SBATCH --time=2:00:00
```

### Adjusting Limits

**GPU Job Limit** (`submit_dist.sh`):
```bash
MAX_GPU_JOBS=48  # Change as needed
```

**Array Size Limit** (`submit_msa_arrays.sh`):
```bash
MAX_ARRAY_SIZE=1000  # SLURM limit
```

**Cycle Duration** (`af3_48hr_cycle.sh`):
```bash
TOTAL_CYCLES=24  # 24 cycles × 2 hours = 48 hours
```

## Troubleshooting

### Common Issues

1. **"Job directory not found"**
   - Verify job exists in `jobs/` or `jobs/human_test_set/`
   - Check job name matches directory exactly
   - Use `get_job_status.sh -v` to see which jobs are missing

2. **MSA jobs not submitting**
   - Run `./monitor_msa_arrays.sh` to check status
   - Verify `.tmp` files exist in base directory
   - Check partition quotas

3. **GPU jobs stuck in queue**
   - Check GPU quota with `squeue -u $USER -p gpu`
   - Verify MSAs are complete: `./get_job_status.sh -stage 2`
   - Run `submit_dist.sh -d` for debug output

4. **Waiting jobs not processing**
   - Ensure MSA arrays have completed
   - Check waiting list: `wc -l waiting_for_msa.csv`
   - Re-run `submit_dist.sh` to process waiting list

5. **No seeds showing for completed jobs**
   - Check output structure: `ls jobs/{job_name}/output/*/seed-*`
   - Verify model.cif files exist in seed directories
   - Use test script if needed: `./test_seed_detection.sh {job_name}`

### Log Files

- **MSA Logs**: `logs/{job_id}_{array_index}_MSA.{out,err}`
- **GPU Logs**: `{job_id}_GPU-only.{out,err}`
- **Cycle Logs**: `af3_cycle_{job_id}_{cycle}.{out,err}`
- **Central Log**: `logged_folding_jobs.csv`

### Recovery Procedures

**After Job Failures**:
```bash
# Re-analyze MSA status
python batch_reuse_msa.py

# Resubmit failed MSA jobs
./submit_msa_arrays.sh msa_array_jobs.csv

# Continue GPU submissions
./submit_dist.sh
```

**After System Interruption**:
```bash
# Check what's running
squeue -u $USER

# Clean up and restart
./cleanup_msa_tmp.sh
./launch_af3.sh
```

## Best Practices

1. **Pre-flight Check**
   - Verify all job directories exist
   - Run `python batch_reuse_msa.py --dry-run` first
   - Check available cluster resources

2. **Regular Monitoring**
   - Check logs every few hours
   - Monitor both partitions for MSA jobs
   - Watch for failed jobs in `logged_folding_jobs.csv`

3. **Resource Management**
   - Don't exceed partition quotas
   - Clean up `.tmp` files after completion
   - Archive completed results regularly

4. **Optimization**
   - Group similar proteins to maximize MSA reuse
   - Balance job distribution between partitions
   - Adjust time limits based on job complexity

## Output Files

### Job Outputs

Each successful job produces:

1. **MSA Stage** (`output_msa/`):
   - `alphafold_input_with_msa.json` - Augmented input
   - `*_data.json` - MSA data files

2. **GPU Stage** (`output/{job_name}/`):
   - `seed-{N}_sample-{0-4}/model.cif` - Structure predictions
   - `confidences.json` - Confidence metrics
   - `summary_confidences.json` - Summary metrics
   - `ranking_scores.csv` - Ranking information

### Pipeline Reports

The monitoring scripts provide:

1. **Pipeline Status** (`pipeline_status.sh`):
   - Current automation cycle
   - Active job counts by type
   - Quick progress estimate

2. **Job Status Report** (`get_job_status.sh`):
   ```
   [STAGE 1] job_name
   [STAGE 2] job_name  
   [STAGE 3] job_name | Seeds: 1,42,64702690
   ```

3. **Detailed Report** (`get_job_status_detailed.sh -e report.csv`):
   ```csv
   job_name,stage,status,completed_seeds,path
   SETD6-CALM1,3,"GPU complete at 2024-11-15 14:30","1,42","/path/to/job"
   ```

### After Pipeline Completion

```bash
# Check final results
tail logged_folding_jobs.csv

# Count successful completions
grep -c ",0$" logged_folding_jobs.csv

# Clean up temporary MSA partition files
./cleanup_msa_tmp.sh

# Sync outputs to organized structure
./sync_organize_outputs.sh

# Upload to Google Drive (requires rclone setup)
./rclone_to_gdrive.sh

# Or use the all-in-one sync script
./sync_all.sh
```

### Output Organization

The sync scripts organize outputs by project:
- `arabidopsis_EEF1A/` - Arabidopsis EEF1A projects
- `ecoli/` - E. coli projects (ecEFTU, ecRPL11)
- `human_EF_KMTs/` - Human elongation factor KMTs
- `yeast/` - Yeast projects (EFM, RKM, SET)
- `SETD6/`, `RPL29/`, `RPL36A/` - Specific protein projects
- `human_test_set/` - Human test set (preserved as-is)

## Quick Reference

### Essential Commands
```bash
# Start pipeline
./launch_af3.sh

# Check status
./pipeline_status.sh              # Quick overview
./get_job_status.sh              # Detailed job stages
./monitor_msa_arrays.sh          # MSA array status

# Debug specific jobs
./get_job_status.sh -f jobs.csv -v
./get_job_status_detailed.sh -stage 2 -p

# Export results
./get_job_status_detailed.sh -e job_report_$(date +%Y%m%d).csv

# Sync outputs
./check_sync_status.sh           # See what's ready
./sync_all.sh                    # Complete sync workflow
./sync_organize_outputs.sh -q    # Local sync (quiet mode)
./rclone_to_gdrive.sh           # Submit Google Drive sync
./check_rclone_status.sh        # Monitor sync progress
```

### Job Stage Reference
- **Stage 1**: Has input, needs MSA → Yellow `[STAGE 1]`
- **Stage 2**: Has MSA, needs GPU → Blue `[STAGE 2]`  
- **Stage 3**: Complete with seeds → Green `[STAGE 3]`

## Support

For issues specific to:
- **Pipeline scripts**: Check this README and script comments
- **AlphaFold3**: Consult AlphaFold3 documentation
- **SLURM**: Contact your HPC support team
- **Email notifications**: Update email in `af3_48hr_cycle.sh`

---

*Pipeline Version: 1.0 | Last Updated: August 2025*
