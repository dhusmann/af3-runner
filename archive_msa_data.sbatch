#!/bin/bash
#SBATCH --job-name=archive_msa
#SBATCH --output=logs/archive_msa_%A_%a.out
#SBATCH --error=logs/archive_msa_%A_%a.err
#SBATCH --partition=hns
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=6:00:00
# Array size will be set dynamically by sync_organize_outputs.sh
# #SBATCH --array=1-50%10

# Archive MSA data files as SLURM array job
# Usage: sbatch archive_msa_data.sbatch
#
# This script automatically discovers job groups and processes them in parallel
# Each array task handles one job group

set -euo pipefail

# Configuration
BASE_DIR="/scratch/groups/ogozani/alphafold3"
SCRIPT_PATH="${BASE_DIR}/archive_msa_data.sh"
JOBS_DIR="${BASE_DIR}/jobs"
LOG_DIR="${BASE_DIR}/logs"

# Create logs directory
mkdir -p "$LOG_DIR"

cd "$BASE_DIR"

# Export environment variables for the archiving script
export COMPRESSION_THREADS=${SLURM_CPUS_PER_TASK:-4}
export DRY_RUN=0

echo "=== MSA Archiving Array Job ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Base directory: $BASE_DIR"
echo "Compression threads: $COMPRESSION_THREADS"
echo

# Function to get all unique job groups
get_job_groups() {
    # Find all job directories and extract unique groups from both locations
    {
        # Main jobs directory
        find "$JOBS_DIR" -maxdepth 1 -type d -name "*-*" | while read job_path; do
            job_name=$(basename "$job_path")
            # Skip human_test_set directory itself
            if [[ "$job_name" != "human_test_set" ]]; then
                # Extract everything before the first hyphen
                echo "${job_name%%-*}"
            fi
        done
        
        # Human test set subdirectory
        if [[ -d "$JOBS_DIR/human_test_set" ]]; then
            find "$JOBS_DIR/human_test_set" -maxdepth 1 -type d -name "*-*" | while read job_path; do
                job_name=$(basename "$job_path")
                # Extract everything before the first hyphen
                echo "${job_name%%-*}"
            done
        fi
    } | sort -u
}

# Get job groups and convert to array
mapfile -t JOB_GROUPS < <(get_job_groups)
TOTAL_GROUPS=${#JOB_GROUPS[@]}

echo "Total job groups discovered: $TOTAL_GROUPS"
echo "Job groups: ${JOB_GROUPS[*]}"
echo

# Check if this array task should process anything
if [[ ${SLURM_ARRAY_TASK_ID} -gt $TOTAL_GROUPS ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} > total groups ($TOTAL_GROUPS)"
    echo "Nothing to process for this task - exiting cleanly"
    exit 0
fi

# Get the job group for this array task (1-indexed)
TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
JOB_GROUP="${JOB_GROUPS[$TASK_INDEX]}"

echo "Processing job group: $JOB_GROUP"
echo "Task index: $TASK_INDEX (array task: ${SLURM_ARRAY_TASK_ID})"
echo

# Check if the archiving script exists
if [[ ! -x "$SCRIPT_PATH" ]]; then
    echo "Error: Archive script not found or not executable: $SCRIPT_PATH"
    exit 1
fi

# Run the archiving script for this specific job group
echo "Executing: $SCRIPT_PATH $JOB_GROUP"
echo "=========================================="
"$SCRIPT_PATH" "$JOB_GROUP"

exit_code=$?

echo
echo "=========================================="
echo "Array task ${SLURM_ARRAY_TASK_ID} completed with exit code: $exit_code"
echo "Job group processed: $JOB_GROUP"
echo "End time: $(date)"

exit $exit_code