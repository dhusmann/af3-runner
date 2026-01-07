# AlphaFold3 Output Sync and Organization

This document describes the **high-performance parallel sync system** for organizing and syncing your AlphaFold3 outputs from the cluster to Google Drive.

## Overview

The system provides a **complete automated workflow** using parallel processing on compute nodes:
1. **Parallel rsync operations**: Organize outputs using 320+ concurrent workers
2. **MSA archiving**: Create deduplicated compressed archives (10-50x compression)
3. **Seed compression**: Parallel compression of seed directories
4. **Google Drive sync**: Upload outputs, scripts, tools, and analysis files
5. **Script backup**: Automatic backup of all your pipeline scripts

**🚀 Performance**: 10-50x faster than the original system, with zero login node impact.

## Key Features

### **Parallel Processing Architecture**
- **320+ concurrent rsync operations** across array tasks
- **All operations on compute nodes** - login nodes remain free
- **Dynamic scaling** - automatically handles any number of jobs
- **Parallel seed compression** across all output directories

### **MSA Management**
- **Smart deduplication** - only archives unique MSA files
- **10-50x compression** through eliminating duplicate MSAs
- **Master index tracking** for easy file location
- **Automated archiving** integrated into workflow
- **Dynamic job group discovery** - processes 486+ job groups automatically
- **Dual structure support** - handles both main jobs and human_test_set directories

### **Complete Backup Solution**
- **Script backup** - all `.sh`, `.py`, `.sbatch`, `.md`, `.conf` files
- **Tools backup** - complete `tools/` directory
- **Analysis backup** - complete `analysis/` directory  
- **Output sync** - organized AlphaFold3 results

### **Automation & Reliability**
- **Zero confirmation prompts** - fully automated execution
- **Google Drive quota handling** with automatic retry
- **Email notifications** when quota limits are hit
- **Job dependency management** between processing steps
- **Incremental syncs** - only new/changed files processed

### **Configuration & Tuning**
- **Configurable parallelism** via `sync_parallel.conf`
- **Resource optimization** for cluster efficiency
- **Dynamic array sizing** based on current job count

## Architecture

### **Processing Flow**
1. **Job Discovery**: Dynamically discovers all jobs needing sync
2. **Parallel Rsync**: Array job organizes outputs using GNU parallel
3. **MSA Archiving**: Array job creates deduplicated MSA archives  
4. **Seed Compression**: Array job compresses seed directories
5. **Google Drive Sync**: Uploads everything including script backups

### **Performance Comparison**
| Operation | Before | After | Improvement |
|-----------|--------|-------|-------------|
| Rsync | 2-4 hours (login node) | 5-15 minutes (compute nodes) | **10-50x faster** |
| MSA Storage | Individual files (slow upload) | Compressed archives | **10-50x reduction** |
| Seed Compression | Sequential | Parallel across directories | **5-20x faster** |
| Overall Workflow | Manual, blocking | Automated, non-blocking | **Complete transformation** |

## Scripts

### **Main Orchestrators**

#### 1. sync_all.sh ⭐
**Complete automated workflow** - the one-command solution:
- Submits all parallel processing jobs
- Submits Google Drive sync with script backup
- No confirmation prompts - fully automated
- Options: `--dry-run`, `--quiet`, `--local-only`

#### 2. sync_organize_outputs.sh
**SLURM job orchestrator** that submits:
- Parallel rsync array job
- MSA archiving array job  
- Seed compression array job
- Options: `--skip-rsync`, `--skip-msa`, `--skip-compress`

#### 3. rclone_to_gdrive.sh
**Enhanced Google Drive sync** that uploads:
- All organized output data
- Complete script backup (`scripts/`)
- Tools directory (`scripts/tools/`)
- Analysis directory (`analysis/`)
- Uses deterministic `--filter` rules for reliable file selection
- 72-hour time limit with quota handling

### **Parallel Processing Workers**

#### 4. sync_organize_rsync.sbatch
**SLURM array job** for parallel rsync operations:
- Dynamic array sizing based on job count
- 16 CPUs, 24GB RAM per task
- GNU parallel for concurrent rsync operations

#### 5. sync_organize_rsync.sh  
**Parallel rsync worker** with organization logic:
- Processes job lists using GNU parallel
- Maintains all original organization patterns
- Progress reporting and error handling

#### 6. archive_msa_data.sbatch
**MSA archiving array job** for parallel processing:
- Groups MSA files by job prefix
- Creates compressed archives with deduplication

#### 7. archive_msa_data.sh
**MSA archiving script** with smart deduplication:
- Finds unique MSA files (skips duplicates) using efficient bash glob patterns
- Creates compressed tar.gz archives by job group
- Maintains master index for file tracking
- Processes 486+ job groups including human_test_set directory
- Handles both flat and nested MSA file structures
- Significant space savings through deduplication
- Fixed hanging issues with optimized file discovery algorithms

#### 8. compress_seeds_array.sbatch
**Parallel seed compression** array job:
- Processes multiple output directories simultaneously
- Creates `seeds.tar.gz` files in each directory
- Much faster than sequential compression

### **Configuration & Discovery**

#### 9. sync_parallel.conf
**Configuration file** for performance tuning:
```bash
JOBS_PER_TASK=500          # Jobs per array task
MAX_ARRAY_TASKS=20         # Maximum concurrent tasks
PARALLEL_WORKERS=16        # Concurrent rsyncs per task
CPUS_PER_TASK=16          # CPU cores per task
MEMORY_PER_TASK="24G"     # RAM per task
```

#### 10. sync_job_discovery.sh
**Dynamic job discovery** system:
- Automatically finds all jobs needing sync
- Calculates optimal array job sizing
- Handles any number of jobs (scales indefinitely)
- Creates job lists for array processing

### **Legacy Scripts (Still Available)**

#### 11. check_sync_status.sh
Status checking utility (less relevant with new architecture)

#### 12. check_rclone_status.sh
Monitor rclone sync jobs and quota status

#### 13. clean_output_dir.sh
Cleanup utility for fresh starts

## Setup

### 1. Make scripts executable
```bash
chmod +x sync_all.sh sync_organize_outputs.sh rclone_to_gdrive.sh
chmod +x sync_organize_rsync.sh sync_job_discovery.sh archive_msa_data.sh
chmod +x *.sbatch
```

### 2. Configure rclone (one-time setup)
```bash
# Install rclone if needed
curl https://rclone.org/install.sh | bash

# Configure Google Drive access
rclone config

# Follow prompts to create remote named 'gozani_labshare_alphafold'
# Test with: rclone lsd gozani_labshare_alphafold:
```

### 3. Optional: Tune performance
Edit `sync_parallel.conf` to adjust parallelism based on your cluster resources and job count.

## Usage

### **Recommended: Complete Automated Workflow**
```bash
# Run everything automatically (recommended)
./sync_all.sh

# Preview all jobs without submitting  
./sync_all.sh --dry-run

# Quiet mode for minimal output
./sync_all.sh --quiet

# Only local processing (skip Google Drive)
./sync_all.sh --local-only
```

### **Individual Components**
```bash
# Submit only local processing jobs
./sync_organize_outputs.sh

# Submit only Google Drive sync
./rclone_to_gdrive.sh

# Test MSA archiving for one group
./archive_msa_data.sh SETD6
```

### **Monitoring Your Jobs**
```bash
# Check all your jobs
squeue -u $USER

# View job details
squeue -j JOB_ID

# Check recent completions
sacct -u $USER --starttime=today
```

## Organization Structure

The system organizes outputs into this structure:

```
output/
├── msa/                   # MSA archives (new)
│   ├── EEF1AKMT1_2024_01_15.tar.gz
│   ├── SETD6_2024_01_15.tar.gz  
│   └── master_index.csv   # File tracking index
├── arabidopsis_EEF1A/     # Jobs with *ateef1a*
├── ecoli/
│   ├── ecEFTU/            # Jobs with *eftu*
│   └── ecRPL11/           # Jobs with *ecrpl*
├── dpEEF1A/               # Jobs with *dpeef1a*
├── RPL29/                 # Jobs with *rpl29*
├── RPL36A/                # Jobs with *rpl36a*
├── SETD6/                 # Jobs with *setd6*
├── human_EF_KMTs/
│   ├── EEF1A/             # Jobs starting with eef1akmt* or mettl13*
│   └── EF2/               # Jobs starting with FAM86*EEF*
├── yeast/
│   ├── EFM/               # Jobs starting with scefm*
│   ├── RKM/               # Jobs starting with scrkm*
│   └── SET/               # Jobs starting with spset*
├── NSD2i/                 # Jobs with aa*[0-9]to[0-9]*
└── human_test_set/        # All human test set jobs
```

## Google Drive Structure

After sync completion, your Google Drive contains:

```
alphafold3/
├── scripts/              # All your pipeline scripts
│   ├── sync_all.sh
│   ├── archive_msa_data.sh
│   ├── [all other scripts]
│   └── tools/            # Contents of tools/ directory
├── analysis/             # Contents of analysis/ directory
│   ├── af3_features.py
│   ├── SETD6.parquet
│   └── [other analysis files]
└── output/               # Organized AlphaFold3 results
    ├── msa/              # MSA archives
    ├── arabidopsis_EEF1A/
    ├── ecoli/
    └── [all other organized results]
```

## MSA Archiving System

### **How It Works**
The MSA archiving system identifies and archives only **unique** MSA files:

1. **Finds original MSAs**: `*_data.json` files in `output_msa/` directories
2. **Skips duplicates**: Ignores `alphafold_input_with_msa.json` (shared MSAs)
3. **Groups by job prefix**: Archives files by job group (e.g., EEF1AKMT1, SETD6)
4. **Creates compressed archives**: `{GROUP}_{DATE}.tar.gz` in `output/msa/`
5. **Maintains master index**: `master_index.csv` for file tracking

### **Space Savings Example**
```
Before: 1000 jobs × 50MB MSA each = 50GB total
After:  10 unique MSAs × 50MB each = 500MB archived
Savings: 99% reduction in storage space
```

### **Query Commands**
```bash
# Find which archive contains a specific job's MSA
grep ",JOB_NAME," output/msa/master_index.csv

# List all unique MSAs in an archive  
grep "^EEF1AKMT1_2024_01_15.tar.gz" output/msa/master_index.csv | cut -d',' -f3

# Count unique MSAs per job group
cut -d',' -f1 output/msa/master_index.csv | cut -d'_' -f1 | sort | uniq -c
```

## What Gets Synced

### **Output Data**
- ✅ All `.cif` model files
- ✅ All `.json` confidence files (except large `*_data.json`)
- ✅ `ranking_scores.csv` files
- ✅ Compressed seed archives (`seeds.tar.gz`)
- ✅ MSA archives in `output/msa/` 
- ❌ Individual `seed-*` directories (now compressed)
- ❌ Large `*_data.json` files (now in MSA archives)
- ❌ `TERMS_OF_USE.md` files

### **Scripts and Analysis** 
- ✅ All script files (`.sh`, `.py`, `.sbatch`, `.md`, `.conf`)
- ✅ Complete `tools/` directory
- ✅ Complete `analysis/` directory
- ❌ Log files (`.out`, `.err`, `.log`, `.tmp`)

## Performance Tuning

### **For Large Job Counts (>20,000 jobs)**
```bash
# Edit sync_parallel.conf
JOBS_PER_TASK=1000
MAX_ARRAY_TASKS=30
PARALLEL_WORKERS=20
```

### **For Limited Cluster Resources**
```bash
# Reduce resource usage
MAX_ARRAY_TASKS=10
PARALLEL_WORKERS=8
CPUS_PER_TASK=8
MEMORY_PER_TASK="16G"
```

### **For Testing with Subsets**
```bash
# Test with limited jobs
MAX_JOBS_LIMIT=1000  # Process first 1000 jobs only
```

## Monitoring

### **Job Status Commands**
```bash
# Check all your running jobs
squeue -u $USER

# Detailed job information
squeue -u $USER -o "%.18i %.12j %.8T %.10M %.6D %R"

# Check specific job
squeue -j JOB_ID

# View completed jobs
sacct -u $USER --starttime=today
```

### **Log Files**
```bash
# Rsync array jobs
ls logs/sync_rsync_*
tail -f logs/sync_rsync_12345_1.out

# MSA archiving
ls logs/archive_msa_*
tail -f logs/archive_msa_12346_1.out

# Seed compression  
ls logs/pack-seeds_*
tail -f logs/pack-seeds_12347_1.out

# Google Drive sync
ls logs/Rclone_*
tail -f logs/Rclone_12348.out
```

### **Progress Tracking**
```bash
# Watch all jobs
watch 'squeue -u $USER'

# Count successful rsyncs
grep "SUCCESS" logs/sync_rsync_*.out | wc -l

# Check MSA archiving progress
grep "Archive created" logs/archive_msa_*.out

# Monitor Google Drive sync steps
grep "STEP [1-4]" logs/Rclone_*.out
```

## Troubleshooting

### **Array Jobs**
```bash
# Check failed array tasks
sacct -j JOB_ID --format=JobID,State,ExitCode

# Resubmit specific failed tasks
sbatch --array=TASK_ID sync_organize_rsync.sbatch

# Check array job limits
sacctmgr show qos
```

### **Performance Issues**
```bash
# Too many jobs for array limits
# Edit sync_parallel.conf and reduce MAX_ARRAY_TASKS

# Jobs taking too long  
# Edit sync_parallel.conf and reduce JOBS_PER_TASK

# Resource conflicts
# Adjust CPUS_PER_TASK and MEMORY_PER_TASK
```

### **Google Drive Issues**
```bash
# Check quota status
grep "quotaExceeded" logs/Rclone_*.out

# Manual retry after quota reset
./rclone_to_gdrive.sh

# Check sync progress by step
grep "=== STEP" logs/Rclone_*.out

# rclone filter warnings (fixed in current version)
# Old versions may show: "Using --filter is recommended instead of both --include and --exclude"
# Solution: Updated to use --filter flags for deterministic rule processing
grep "filter is recommended" logs/Rclone_*.out  # Should show no results
```

### **MSA Archiving Issues**
```bash
# Check master index
head -5 output/msa/master_index.csv

# Verify archive integrity
tar -tzf output/msa/EEF1AKMT1_*.tar.gz | head

# Find jobs without unique MSAs
grep "No unique MSA found" logs/archive_msa_*.out

# Script hanging during MSA archiving
# Solution: The script uses optimized bash glob patterns instead of nested for loops
# If hanging occurs, check that archive_msa_data.sh uses direct directory iteration

# Check if all job groups are being processed
grep "Processing job group" logs/archive_msa_*.out | wc -l  # Should show 486+ groups

# Verify both main and human_test_set directories are included
ls -la jobs/ | grep human_test_set  # Should show human_test_set directory
```

### **Recovery Commands**
```bash
# Cancel all running jobs
scancel -u $USER

# Clean up and restart
./clean_output_dir.sh  # if needed
./sync_all.sh

# Check system status
df -h /scratch/groups/ogozani/alphafold3/
```

## Advanced Usage

### **Partial Runs**
```bash
# Skip specific components
./sync_organize_outputs.sh --skip-msa      # Skip MSA archiving
./sync_organize_outputs.sh --skip-compress # Skip seed compression  
./sync_organize_outputs.sh --rsync-only    # Only rsync

# Process specific MSA groups
./archive_msa_data.sh EEF1AKMT1  # Single group
export MAX_JOBS_LIMIT=1000       # Limit for testing
```

### **Configuration Examples**
```bash
# High-performance cluster setup
JOBS_PER_TASK=2000
MAX_ARRAY_TASKS=50
PARALLEL_WORKERS=32
CPUS_PER_TASK=32
MEMORY_PER_TASK="64G"

# Conservative resource usage
JOBS_PER_TASK=200
MAX_ARRAY_TASKS=5  
PARALLEL_WORKERS=4
CPUS_PER_TASK=4
MEMORY_PER_TASK="8G"
```

## Storage Estimates

### **Space Savings**
- **MSA deduplication**: 10-50x reduction in MSA storage
- **Seed compression**: ~50-70% reduction in seed storage  
- **Overall improvement**: 30-80% less storage needed
- **Upload speed**: 5-10x faster with compressed archives

### **File Size References**
- Model files (`.cif`): ~1-5 MB each
- Confidence files (`.json`): ~100 KB each  
- Seed archive (`seeds.tar.gz`): ~10-50 MB per job
- MSA archive: ~10-500 MB per unique group
- Script backup: ~1-10 MB total

## Migration from Old System

If you were using the previous sync system:

1. **Backup existing configs** (automatic backups created)
2. **Run new system**: `./sync_all.sh --dry-run` first
3. **Monitor logs** for any issues
4. **Enjoy the performance boost**! 🚀

The new system is fully backward compatible and maintains the same organization structure.

## Next Steps

After your sync completes, you can:

1. **Monitor progress** with the provided monitoring commands
2. **Share specific folders** from Google Drive with collaborators  
3. **Download results** to local machines for analysis
4. **Use Google Colab** for structure visualization
5. **Set up automated workflows** leveraging the script backups
6. **Scale to larger datasets** - the system handles growth automatically

---

**🎉 Congratulations!** You now have a high-performance, fully automated sync system that scales with your research and keeps your login nodes free for interactive work.