#!/bin/bash

# sync_job_discovery.sh - Dynamic job discovery and chunking for parallel sync
#
# This script discovers all jobs that need syncing and creates job lists
# for array processing. It's designed to handle any number of jobs dynamically.
#
# Usage:
#   ./sync_job_discovery.sh [action] [options]
#
# Actions:
#   count     - Count total jobs that need syncing
#   list      - List all jobs that need syncing
#   chunks    - Create job chunks for array processing
#   info      - Show discovery information

set -euo pipefail

# Load configuration
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
CONFIG_FILE="${SCRIPT_DIR}/sync_parallel.conf"

if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Warning: Configuration file not found: $CONFIG_FILE"
    echo "Using default settings"
    
    # Default settings
    BASE_DIR="/scratch/groups/ogozani/alphafold3"
    JOBS_DIR="${BASE_DIR}/jobs"
    OUTPUT_DIR="${BASE_DIR}/output"
    JOBS_PER_TASK=500
    MAX_JOBS_LIMIT=0
    DEBUG_MODE=false
fi

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Action from command line
ACTION="${1:-info}"

# Function to discover all jobs that need syncing (optimized for large datasets)
discover_jobs() {
    local jobs_list=$(mktemp)
    
    if [[ $DEBUG_MODE == "true" ]]; then
        echo "Discovering jobs in: $JOBS_DIR" >&2
    fi
    
    # More efficient approach - use find with -exec instead of nested while loops
    # Process top-level jobs (excluding human_test_set)
    find "$JOBS_DIR" -maxdepth 2 -type d -path "*/output" ! -path "*/human_test_set/*" -exec bash -c '
        output_dir="$1"
        job_dir="$(dirname "$output_dir")"
        for subdir in "$output_dir"/*/; do
            if [[ -d "$subdir" ]]; then
                echo "TOPLEVEL|${job_dir}|${subdir%/}"
            fi
        done
    ' _ {} \; >> "$jobs_list"
    
    # Process human_test_set jobs
    if [[ -d "${JOBS_DIR}/human_test_set" ]]; then
        find "${JOBS_DIR}/human_test_set" -maxdepth 2 -type d -path "*/output" -exec bash -c '
            output_dir="$1"
            job_dir="$(dirname "$output_dir")"
            for subdir in "$output_dir"/*/; do
                if [[ -d "$subdir" ]]; then
                    echo "HUMAN_TEST_SET|${job_dir}|${subdir%/}"
                fi
            done
        ' _ {} \; >> "$jobs_list"
    fi
    
    # Apply job limit if set (for testing)
    if [[ $MAX_JOBS_LIMIT -gt 0 ]]; then
        head -n "$MAX_JOBS_LIMIT" "$jobs_list" > "${jobs_list}.limited"
        mv "${jobs_list}.limited" "$jobs_list"
    fi
    
    echo "$jobs_list"
}

# Function to count jobs
count_jobs() {
    local jobs_list=$(discover_jobs)
    local count=$(wc -l < "$jobs_list")
    rm -f "$jobs_list"
    echo "$count"
}

# Function to list jobs
list_jobs() {
    local jobs_list=$(discover_jobs)
    cat "$jobs_list"
    rm -f "$jobs_list"
}

# Function to get job chunks for array processing
get_job_chunks() {
    local jobs_list=$(discover_jobs)
    local total_jobs=$(wc -l < "$jobs_list")
    
    if [[ $total_jobs -eq 0 ]]; then
        echo "No jobs found to process" >&2
        rm -f "$jobs_list"
        return 1
    fi
    
    # Calculate array size
    local array_size=$(( (total_jobs + JOBS_PER_TASK - 1) / JOBS_PER_TASK ))
    
    # Ensure we don't exceed max array tasks
    if [[ $array_size -gt $MAX_ARRAY_TASKS ]]; then
        echo "Warning: Calculated array size ($array_size) exceeds maximum ($MAX_ARRAY_TASKS)" >&2
        echo "Consider increasing JOBS_PER_TASK or MAX_ARRAY_TASKS in config" >&2
        array_size=$MAX_ARRAY_TASKS
        JOBS_PER_TASK=$(( (total_jobs + array_size - 1) / array_size ))
    fi
    
    echo "TOTAL_JOBS=$total_jobs"
    echo "ARRAY_SIZE=$array_size"
    echo "JOBS_PER_TASK=$JOBS_PER_TASK"
    echo "JOBS_LIST=$jobs_list"
    
    if [[ $DEBUG_MODE == "true" ]]; then
        echo "Job distribution:" >&2
        for ((i=1; i<=array_size; i++)); do
            local start_line=$(( (i - 1) * JOBS_PER_TASK + 1 ))
            local end_line=$(( i * JOBS_PER_TASK ))
            [[ $end_line -gt $total_jobs ]] && end_line=$total_jobs
            local chunk_size=$(( end_line - start_line + 1 ))
            echo "  Task $i: lines $start_line-$end_line ($chunk_size jobs)" >&2
        done
    fi
}

# Function to get jobs for specific array task
get_array_task_jobs() {
    local array_task_id="$1"
    local jobs_list="$2"
    local total_jobs=$(wc -l < "$jobs_list")
    local array_size=$(( (total_jobs + JOBS_PER_TASK - 1) / JOBS_PER_TASK ))
    
    # Calculate line range for this task
    local start_line=$(( (array_task_id - 1) * JOBS_PER_TASK + 1 ))
    local end_line=$(( array_task_id * JOBS_PER_TASK ))
    [[ $end_line -gt $total_jobs ]] && end_line=$total_jobs
    
    if [[ $start_line -le $total_jobs ]]; then
        sed -n "${start_line},${end_line}p" "$jobs_list"
    fi
}

# Function to show discovery info
show_info() {
    echo -e "${BLUE}AlphaFold3 Job Discovery Information${NC}"
    echo "=================================="
    echo "Base directory: $BASE_DIR"
    echo "Jobs directory: $JOBS_DIR"
    echo "Output directory: $OUTPUT_DIR"
    echo
    
    local jobs_list=$(discover_jobs)
    local total_jobs=$(wc -l < "$jobs_list")
    local array_size=$(( (total_jobs + JOBS_PER_TASK - 1) / JOBS_PER_TASK ))
    
    # Count job types
    local toplevel_jobs=$(grep -c "^TOPLEVEL" "$jobs_list" || echo 0)
    local human_test_jobs=$(grep -c "^HUMAN_TEST_SET" "$jobs_list" || echo 0)
    
    echo "Job Statistics:"
    echo "  Total jobs: $total_jobs"
    echo "  Top-level jobs: $toplevel_jobs"
    echo "  Human test set jobs: $human_test_jobs"
    echo
    
    if [[ $MAX_JOBS_LIMIT -gt 0 ]]; then
        echo -e "${YELLOW}Note: Job limit active - showing first $MAX_JOBS_LIMIT jobs only${NC}"
        echo
    fi
    
    echo "Array Configuration:"
    echo "  Jobs per task: $JOBS_PER_TASK"
    echo "  Calculated array size: $array_size"
    echo "  Max array tasks: $MAX_ARRAY_TASKS"
    echo "  Parallel workers: $PARALLEL_WORKERS"
    
    if [[ $array_size -gt $MAX_ARRAY_TASKS ]]; then
        echo -e "${YELLOW}  Warning: Array size exceeds maximum${NC}"
        local adjusted_jobs_per_task=$(( (total_jobs + MAX_ARRAY_TASKS - 1) / MAX_ARRAY_TASKS ))
        echo "  Suggested JOBS_PER_TASK: $adjusted_jobs_per_task"
    fi
    
    echo
    echo "Estimated Resources:"
    echo "  Total CPU cores: $(( array_size * PARALLEL_WORKERS ))"
    echo "  Total memory: $(( array_size )) × $MEMORY_PER_TASK"
    
    rm -f "$jobs_list"
}

# Main execution
case "$ACTION" in
    "count")
        count_jobs
        ;;
    "list")
        list_jobs
        ;;
    "chunks")
        get_job_chunks
        ;;
    "task")
        # Get jobs for specific array task
        # Usage: sync_job_discovery.sh task TASK_ID JOBS_LIST_FILE
        if [[ $# -lt 3 ]]; then
            echo "Usage: $0 task TASK_ID JOBS_LIST_FILE" >&2
            exit 1
        fi
        get_array_task_jobs "$2" "$3"
        ;;
    "info")
        show_info
        ;;
    *)
        echo "Unknown action: $ACTION" >&2
        echo "Valid actions: count, list, chunks, task, info" >&2
        exit 1
        ;;
esac