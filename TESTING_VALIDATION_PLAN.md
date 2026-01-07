# AlphaFold3 Output Synchronization - Testing & Validation Plan

## Overview

This document outlines the comprehensive testing plan for the improved AlphaFold3 output synchronization workflow, including MSA archiving with deduplication and parallel processing improvements.

## Pre-Testing Setup

### 1. Environment Verification
```bash
# Verify all scripts are executable
ls -la /scratch/groups/ogozani/alphafold3/archive_msa_data.sh
ls -la /scratch/groups/ogozani/alphafold3/archive_msa_data.sbatch
ls -la /scratch/groups/ogozani/alphafold3/compress_seeds_array.sbatch
ls -la /scratch/groups/ogozani/alphafold3/sync_organize_outputs.sh
ls -la /scratch/groups/ogozani/alphafold3/rclone_to_gdrive.sh

# Verify required tools are available
which pigz
which rclone
which bc
```

### 2. Backup Existing Configuration
```bash
# Create backup of original scripts
mkdir -p /scratch/groups/ogozani/alphafold3/backup_$(date +%Y%m%d)
cp compress_seeds.sbatch backup_*/
cp sync_organize_outputs.sh backup_*/sync_organize_outputs.sh.orig
cp rclone_to_gdrive.sh backup_*/rclone_to_gdrive.sh.orig
```

## Testing Phase 1: MSA Archiving Script Validation

### Test 1.1: Dry Run MSA Discovery
```bash
cd /scratch/groups/ogozani/alphafold3
export DRY_RUN=1
./archive_msa_data.sh
```
**Expected Results:**
- Script discovers all job groups
- Shows count of unique MSA files vs shared MSAs
- Shows what archives would be created
- No files actually created

### Test 1.2: Single Job Group Archive
```bash
# Test with a small job group first (e.g., SETD6)
export DRY_RUN=0
./archive_msa_data.sh SETD6
```
**Expected Results:**
- Creates archive in `/scratch/groups/ogozani/alphafold3/output/msa/`
- Archive named `SETD6_YYYY_MM_DD.tar.gz`
- Master index updated with entries
- Original MSA files remain untouched

### Test 1.3: Archive Verification
```bash
# Verify archive contents
cd /scratch/groups/ogozani/alphafold3/output/msa
tar -tzf SETD6_*.tar.gz | head -10
tar -tzf SETD6_*.tar.gz | wc -l

# Verify master index
head -5 master_index.csv
grep "SETD6" master_index.csv | wc -l
```

### Test 1.4: Deduplication Verification
```bash
# Test query commands
grep ",SETD6-" master_index.csv
grep "^SETD6_" master_index.csv | cut -d',' -f3 | sort -u
```

### Test 1.5: Re-run Protection
```bash
# Run again to verify already-archived files are skipped
./archive_msa_data.sh SETD6
```
**Expected Results:**
- Script reports files already archived
- No duplicate entries in master index
- No duplicate archives created

## Testing Phase 2: SLURM Array Job Testing

### Test 2.1: MSA Array Job - Small Scale
```bash
# Submit with limited array size for testing
sbatch --array=1-3%2 archive_msa_data.sbatch
```
**Expected Results:**
- 3 array tasks created
- Only 2 running simultaneously (%2 limit)
- Each task processes one job group
- Logs created in `logs/archive_msa_*.out`

### Test 2.2: Seeds Compression Array Job
```bash
# Test the new array-based compression
sbatch --array=1-5%3 compress_seeds_array.sbatch
```
**Expected Results:**
- Processes 5 directories in parallel (max 3 at once)
- Creates `seeds.tar.gz` files in each output directory
- Proper resource utilization (8 CPUs, 16GB per task)

### Test 2.3: Job Monitoring
```bash
# Monitor jobs
squeue -u $USER
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed

# Check logs
tail -f logs/archive_msa_*_*.out
tail -f logs/pack-seeds_*_*.out
```

## Testing Phase 3: Integrated Workflow Testing

### Test 3.1: Updated sync_organize_outputs.sh
```bash
# Test with dry run first
./sync_organize_outputs.sh --dry-run
```
**Expected Results:**
- Shows what would be synced
- Reports MSA and compression jobs would be submitted
- No actual job submission in dry-run mode

### Test 3.2: Full Workflow Test (Small Scale)
```bash
# Run actual sync with job submission
./sync_organize_outputs.sh --quiet
```
**Expected Results:**
- Jobs synced as before
- MSA archiving job submitted
- Seeds compression array job submitted
- Job IDs reported
- `output/msa/` directory created

### Test 3.3: Verify Job Dependencies
```bash
# Check that both jobs are running
squeue -u $USER -o "%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R"

# Verify resource allocation
sstat -j JOB_ID --format=JobID,AveCPU,AveRSS,AveVMSize
```

## Testing Phase 4: rclone Integration Testing

### Test 4.1: rclone Configuration Verification
```bash
# Test rclone connectivity (dry run)
rclone lsd gozani_labshare_alphafold:alphafold3/output/ --dry-run
```

### Test 4.2: Upload Analysis
```bash
# Run with updated file analysis
./rclone_to_gdrive.sh
```
**Expected Results:**
- Reports seed archives count
- Reports MSA archives count and size
- Includes MSA directory in sync
- Excludes uncompressed seed directories

### Test 4.3: Selective Upload Test
```bash
# Test that only correct files are uploaded
# Create small test structure
mkdir -p test_output/msa
echo "test" > test_output/test.txt
echo "test_msa" > test_output/msa/test_msa.tar.gz
mkdir -p test_output/job/seed-1
echo "should_not_sync" > test_output/job/seed-1/file.txt
echo "should_sync" > test_output/job/seeds.tar.gz

# Test rclone command manually
rclone copy -v test_output/ gdrive_remote:test/ \
    --exclude "**/seed-*/**" \
    --exclude "**/seed-*" \
    --include "**/seeds.tar.gz" \
    --include "**/seeds.tar" \
    --include "msa/**" \
    --dry-run
```

## Testing Phase 5: Performance & Scale Testing

### Test 5.1: Full MSA Archive Creation
```bash
# Run full MSA archiving for all job groups
sbatch archive_msa_data.sbatch
```

### Test 5.2: Compression Ratio Analysis
```bash
# Analyze compression effectiveness
cd /scratch/groups/ogozani/alphafold3/output/msa
for archive in *.tar.gz; do
    size=$(stat -c%s "$archive")
    echo "$archive: $(numfmt --to=iec $size)"
done | sort -k2 -hr
```

### Test 5.3: Parallel Processing Performance
```bash
# Compare old vs new compression times
time sbatch --wait compress_seeds.sbatch  # Original
time sbatch --wait compress_seeds_array.sbatch  # New array job
```

### Test 5.4: Storage Space Verification
```bash
# Calculate total space savings
original_msa_size=$(find /scratch/groups/ogozani/alphafold3/jobs -name "*_data.json" -exec stat -c%s {} \; | awk '{sum+=$1} END {print sum}')
archived_msa_size=$(find /scratch/groups/ogozani/alphafold3/output/msa -name "*.tar.gz" -exec stat -c%s {} \; | awk '{sum+=$1} END {print sum}')

echo "Original MSA size: $(numfmt --to=iec $original_msa_size)"
echo "Archived MSA size: $(numfmt --to=iec $archived_msa_size)"
echo "Space saved: $(numfmt --to=iec $((original_msa_size - archived_msa_size)))"
```

## Testing Phase 6: Error Handling & Recovery

### Test 6.1: Partial Archive Recovery
```bash
# Simulate interrupted archiving
# Create incomplete master index and test recovery
echo "archive1.tar.gz,job1,file1.json,1000,2024-01-01" >> output/msa/master_index.csv
./archive_msa_data.sh
```

### Test 6.2: Storage Full Simulation
```bash
# Test behavior when output directory is full
# (Requires careful setup in test environment)
```

### Test 6.3: Network Interruption Recovery
```bash
# Test rclone resumability after interruption
# (Monitor and interrupt rclone job manually)
```

## Validation Criteria

### Success Criteria
1. **MSA Archiving:**
   - [ ] All unique MSA files identified and archived
   - [ ] No duplicate MSAs in archives (deduplication working)
   - [ ] Master index accurately tracks all archived files
   - [ ] Archives are valid and extractable
   - [ ] Space reduction of 10-50x achieved

2. **Parallel Processing:**
   - [ ] Array jobs distribute work correctly
   - [ ] Resource utilization improved over serial processing
   - [ ] No race conditions or file conflicts
   - [ ] All jobs complete successfully

3. **Integration:**
   - [ ] sync_organize_outputs.sh submits both job types
   - [ ] rclone syncs MSA archives to Google Drive
   - [ ] No disruption to existing workflow
   - [ ] Proper error reporting and logging

4. **Data Integrity:**
   - [ ] No original files modified or deleted
   - [ ] Archive contents match source files exactly
   - [ ] Master index allows accurate file location
   - [ ] All compressed files are valid

### Performance Benchmarks
- MSA archiving should complete within 6 hours
- Seed compression should be faster than serial version
- Total storage reduction of at least 30% for MSA files
- Upload speed improvement of at least 2x for MSA data

## Post-Testing Cleanup

### Successful Test Completion
```bash
# Update documentation
# Create final README update
# Archive test logs
mkdir -p logs/testing_$(date +%Y%m%d)
mv logs/archive_msa_* logs/testing_*/
mv logs/pack-seeds_* logs/testing_*/
```

### Failed Test Recovery
```bash
# Restore backups if needed
cp backup_*/sync_organize_outputs.sh.orig sync_organize_outputs.sh
cp backup_*/rclone_to_gdrive.sh.orig rclone_to_gdrive.sh
cp backup_*/compress_seeds.sbatch .

# Clean up test artifacts
rm -rf output/msa/  # If needed
```

## Long-term Monitoring

### Weekly Checks
- Verify master index integrity
- Check archive file sizes and compression ratios
- Monitor SLURM job success rates
- Validate Google Drive sync completeness

### Monthly Reviews
- Analyze storage space trends
- Review job performance metrics
- Update array job size limits if needed
- Clean up old log files

### Troubleshooting Reference
Common issues and solutions:
1. **Array jobs exceed available groups:** Adjust `--array=1-N` in sbatch files
2. **pigz not available:** Scripts fall back to gzip automatically
3. **rclone quota exceeded:** Built-in retry mechanism handles this
4. **Master index corruption:** Rebuild from archive contents
5. **SLURM resource limits:** Adjust memory/CPU requirements in sbatch files