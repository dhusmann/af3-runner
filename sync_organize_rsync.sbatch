#!/bin/bash
#SBATCH --job-name=sync_rsync_array
#SBATCH --output=logs/sync_rsync_%A_%a.out
#SBATCH --error=logs/sync_rsync_%A_%a.err
#SBATCH --partition=hns
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=4:00:00

# Dynamic SLURM array job for parallel rsync operations
# This script is submitted by sync_organize_outputs.sh with dynamically calculated array size
# Each array task processes a chunk of jobs using GNU parallel

set -euo pipefail

# Configuration
BASE_DIR="/scratch/groups/ogozani/alphafold3"
SCRIPT_DIR="$BASE_DIR"
RSYNC_WORKER="${SCRIPT_DIR}/sync_organize_rsync.sh"
JOB_DISCOVERY="${SCRIPT_DIR}/sync_job_discovery.sh"
CONFIG_FILE="${SCRIPT_DIR}/sync_parallel.conf"
LOG_DIR="${BASE_DIR}/logs"

# Ensure we're in the base directory
cd "$BASE_DIR"

# Load configuration
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Error: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

# Create logs directory
mkdir -p "$LOG_DIR"

echo "=== AlphaFold3 Parallel Rsync Array Job ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "Base directory: $BASE_DIR"
echo "CPUs allocated: ${SLURM_CPUS_PER_TASK}"
echo "Memory allocated: $MEMORY_PER_TASK"
echo

# Verify required scripts exist
if [[ ! -x "$RSYNC_WORKER" ]]; then
    echo "Error: Rsync worker script not found or not executable: $RSYNC_WORKER"
    exit 1
fi

if [[ ! -x "$JOB_DISCOVERY" ]]; then
    echo "Error: Job discovery script not found or not executable: $JOB_DISCOVERY"
    exit 1
fi

# Check for required environment variables passed from orchestrator
if [[ -z "${SYNC_JOBS_LIST:-}" ]]; then
    echo "Error: SYNC_JOBS_LIST environment variable not set"
    echo "This script should be called by sync_organize_outputs.sh"
    exit 1
fi

if [[ -z "${SYNC_TOTAL_JOBS:-}" ]]; then
    echo "Error: SYNC_TOTAL_JOBS environment variable not set"
    exit 1
fi

if [[ -z "${SYNC_JOBS_PER_TASK:-}" ]]; then
    echo "Error: SYNC_JOBS_PER_TASK environment variable not set"
    exit 1
fi

# Verify jobs list file exists
if [[ ! -f "$SYNC_JOBS_LIST" ]]; then
    echo "Error: Jobs list file not found: $SYNC_JOBS_LIST"
    exit 1
fi

echo "Jobs list file: $SYNC_JOBS_LIST"
echo "Total jobs in cluster: $SYNC_TOTAL_JOBS"
echo "Jobs per array task: $SYNC_JOBS_PER_TASK"
echo "Parallel workers per task: $PARALLEL_WORKERS"
echo

# Calculate job range for this array task
START_LINE=$(( (SLURM_ARRAY_TASK_ID - 1) * SYNC_JOBS_PER_TASK + 1 ))
END_LINE=$(( SLURM_ARRAY_TASK_ID * SYNC_JOBS_PER_TASK ))

# Don't exceed total jobs
if [[ $END_LINE -gt $SYNC_TOTAL_JOBS ]]; then
    END_LINE=$SYNC_TOTAL_JOBS
fi

# Check if this array task should process anything
if [[ $START_LINE -gt $SYNC_TOTAL_JOBS ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} > total jobs ($SYNC_TOTAL_JOBS)"
    echo "Nothing to process for this task - exiting cleanly"
    exit 0
fi

CHUNK_SIZE=$(( END_LINE - START_LINE + 1 ))

echo "Processing job range: $START_LINE to $END_LINE ($CHUNK_SIZE jobs)"
echo "Effective parallelism: $PARALLEL_WORKERS concurrent rsyncs"
echo

# Extract jobs for this array task
TASK_JOBS_FILE=$(mktemp)
sed -n "${START_LINE},${END_LINE}p" "$SYNC_JOBS_LIST" > "$TASK_JOBS_FILE"

if [[ ! -s "$TASK_JOBS_FILE" ]]; then
    echo "No jobs found in range $START_LINE-$END_LINE"
    rm -f "$TASK_JOBS_FILE"
    exit 0
fi

echo "Jobs for this task:"
head -5 "$TASK_JOBS_FILE"
if [[ $CHUNK_SIZE -gt 5 ]]; then
    echo "... and $(( CHUNK_SIZE - 5 )) more jobs"
fi
echo

# Run the rsync worker with parallel processing
echo "Starting parallel rsync processing..."
echo "=========================================="

# Set environment for worker script
export PARALLEL_WORKERS
export DEBUG_MODE

# Execute the rsync worker
if "$RSYNC_WORKER" process_file "$TASK_JOBS_FILE" "$PARALLEL_WORKERS"; then
    exit_code=0
    echo "=========================================="
    echo "Array task ${SLURM_ARRAY_TASK_ID} completed successfully"
else
    exit_code=$?
    echo "=========================================="
    echo "Array task ${SLURM_ARRAY_TASK_ID} completed with errors"
fi

# Cleanup
rm -f "$TASK_JOBS_FILE"

echo "Job range processed: $START_LINE-$END_LINE"
echo "Chunk size: $CHUNK_SIZE jobs"
echo "End time: $(date)"
echo "Exit code: $exit_code"

# Report resource usage
if command -v sacct >/dev/null 2>&1; then
    echo
    echo "Resource usage:"
    sacct -j "${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}" --format=JobID,MaxRSS,MaxVMSize,AveCPU,Elapsed 2>/dev/null || true
fi

exit $exit_code