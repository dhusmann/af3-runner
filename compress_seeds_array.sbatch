#!/bin/bash
#SBATCH --job-name=pack-seeds-array
#SBATCH --output=logs/pack-seeds_%A_%a.out
#SBATCH --error=logs/pack-seeds_%A_%a.err
#SBATCH --partition=hns
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-20%8

# Compress seeds as SLURM array job for parallel processing
# Usage: sbatch compress_seeds_array.sbatch
#
# This script automatically discovers output directories and processes them in parallel
# Each array task handles one output directory

set -euo pipefail

# Configuration
BASE_DIR="/scratch/groups/ogozani/alphafold3"
OUTPUT_DIR="${BASE_DIR}/output"
SCRIPT_PATH="${BASE_DIR}/compress_seeds.sh"
LOG_DIR="${BASE_DIR}/logs"

# Create logs directory
mkdir -p "$LOG_DIR"

cd "$BASE_DIR"

# Export environment variables for compression script
export THREADS=${SLURM_CPUS_PER_TASK:-8}
export COMPRESSION_LEVEL=7
export PARALLEL_DIRS=1
export COMPRESSOR=auto

echo "=== Seed Compression Array Job ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Base directory: $BASE_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Compression threads: $THREADS"
echo

# Function to get all output directories that need compression
get_output_directories() {
    local dirs=()
    
    # Add all organized project directories
    for dir in \
        "${OUTPUT_DIR}/arabidopsis_EEF1A" \
        "${OUTPUT_DIR}/ecoli/ecEFTU" \
        "${OUTPUT_DIR}/ecoli/ecRPL11" \
        "${OUTPUT_DIR}/dpEEF1A" \
        "${OUTPUT_DIR}/RPL29" \
        "${OUTPUT_DIR}/RPL36A" \
        "${OUTPUT_DIR}/SETD6" \
        "${OUTPUT_DIR}/human_EF_KMTs/EEF1A" \
        "${OUTPUT_DIR}/human_EF_KMTs/EF2" \
        "${OUTPUT_DIR}/yeast/EFM" \
        "${OUTPUT_DIR}/yeast/RKM" \
        "${OUTPUT_DIR}/yeast/SET" \
        "${OUTPUT_DIR}/NSD2i" \
        "${OUTPUT_DIR}/human_test_set"; do
        
        if [[ -d "$dir" ]]; then
            dirs+=("$dir")
        fi
    done
    
    # Add any other top-level directories (unmatched jobs)
    for dir in "${OUTPUT_DIR}"*/; do
        if [[ -d "$dir" ]]; then
            dir_name=$(basename "$dir")
            # Skip if already included in organized directories
            case "$dir_name" in
                arabidopsis_EEF1A|ecoli|dpEEF1A|RPL29|RPL36A|SETD6|human_EF_KMTs|yeast|NSD2i|human_test_set|msa)
                    ;;
                *)
                    dirs+=("$dir")
                    ;;
            esac
        fi
    done
    
    printf '%s\n' "${dirs[@]}"
}

# Get directories to process and convert to array
mapfile -t DIRECTORIES < <(get_output_directories)
TOTAL_DIRS=${#DIRECTORIES[@]}

echo "Total directories discovered: $TOTAL_DIRS"
echo "Directories to process:"
for dir in "${DIRECTORIES[@]}"; do
    echo "  - $dir"
done
echo

# Check if this array task should process anything
if [[ ${SLURM_ARRAY_TASK_ID} -gt $TOTAL_DIRS ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} > total directories ($TOTAL_DIRS)"
    echo "Nothing to process for this task - exiting cleanly"
    exit 0
fi

# Get the directory for this array task (1-indexed)
TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
TARGET_DIR="${DIRECTORIES[$TASK_INDEX]}"

echo "Processing directory: $TARGET_DIR"
echo "Task index: $TASK_INDEX (array task: ${SLURM_ARRAY_TASK_ID})"
echo

# Check if the compression script exists
if [[ ! -x "$SCRIPT_PATH" ]]; then
    echo "Error: Compression script not found or not executable: $SCRIPT_PATH"
    exit 1
fi

# Check if target directory exists
if [[ ! -d "$TARGET_DIR" ]]; then
    echo "Error: Target directory does not exist: $TARGET_DIR"
    exit 1
fi

# Count subdirectories with seed-* folders before compression
BEFORE_COUNT=$(find "$TARGET_DIR" -name "seed-*" -type d | wc -l)
echo "Found $BEFORE_COUNT seed directories to potentially compress"

# Run the compression script for this specific directory
echo "Executing: $SCRIPT_PATH $TARGET_DIR"
echo "=========================================="
"$SCRIPT_PATH" "$TARGET_DIR"

exit_code=$?

# Count archives created
AFTER_COUNT=$(find "$TARGET_DIR" -name "seeds.tar.gz" -o -name "seeds.tar" | wc -l)

echo
echo "=========================================="
echo "Array task ${SLURM_ARRAY_TASK_ID} completed with exit code: $exit_code"
echo "Directory processed: $TARGET_DIR"
echo "Archives created/verified: $AFTER_COUNT"
echo "End time: $(date)"

exit $exit_code