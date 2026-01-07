#!/bin/bash

# sync_organize_outputs.sh - AlphaFold3 Output Sync Orchestrator
# 
# This script orchestrates all sync operations by submitting SLURM jobs:
# 1. Parallel rsync for local organization (compute nodes)
# 2. MSA archiving (parallel array job)
# 3. Seed compression (parallel array job)
#
# All heavy I/O operations now run on compute nodes instead of login nodes.

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Base directories
BASE_DIR="/scratch/groups/ogozani/alphafold3"
JOBS_DIR="${BASE_DIR}/jobs"
OUTPUT_DIR="${BASE_DIR}/output"
SCRIPT_DIR="$BASE_DIR"

# Required scripts
JOB_DISCOVERY="${SCRIPT_DIR}/sync_job_discovery.sh"
RSYNC_SBATCH="${SCRIPT_DIR}/sync_organize_rsync.sbatch"
MSA_SBATCH="${SCRIPT_DIR}/archive_msa_data.sbatch"
COMPRESS_SBATCH="${SCRIPT_DIR}/compress_seeds_array.sbatch"
CONFIG_FILE="${SCRIPT_DIR}/sync_parallel.conf"

# Load configuration
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo -e "${RED}Error: Configuration file not found: $CONFIG_FILE${NC}"
    exit 1
fi

# Default options
DRY_RUN=0
QUIET=0
VERBOSE=0
SKIP_RSYNC=0
SKIP_MSA=0
SKIP_COMPRESS=0

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--dry-run)
            DRY_RUN=1
            shift
            ;;
        -q|--quiet)
            QUIET=1
            shift
            ;;
        -v|--verbose)
            VERBOSE=1
            QUIET=0
            shift
            ;;
        --skip-rsync)
            SKIP_RSYNC=1
            shift
            ;;
        --skip-msa)
            SKIP_MSA=1
            shift
            ;;
        --skip-compress)
            SKIP_COMPRESS=1
            shift
            ;;
        --rsync-only)
            SKIP_MSA=1
            SKIP_COMPRESS=1
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  -n, --dry-run       Preview job submissions without executing"
            echo "  -q, --quiet         Minimal output (summary only)"
            echo "  -v, --verbose       Detailed output"
            echo "  --skip-rsync        Skip rsync job submission"
            echo "  --skip-msa          Skip MSA archiving job submission"
            echo "  --skip-compress     Skip seed compression job submission"
            echo "  --rsync-only        Only submit rsync job"
            echo "  -h, --help          Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Verify required scripts exist
check_dependencies() {
    local missing_deps=0
    
    for script in "$JOB_DISCOVERY" "$RSYNC_SBATCH"; do
        if [[ ! -x "$script" ]]; then
            echo -e "${RED}Error: Required script not found or not executable: $script${NC}"
            missing_deps=1
        fi
    done
    
    if [[ $SKIP_MSA -eq 0 && ! -x "$MSA_SBATCH" ]]; then
        echo -e "${RED}Error: MSA archiving script not found: $MSA_SBATCH${NC}"
        missing_deps=1
    fi
    
    if [[ $SKIP_COMPRESS -eq 0 && ! -x "$COMPRESS_SBATCH" ]]; then
        echo -e "${RED}Error: Seed compression script not found: $COMPRESS_SBATCH${NC}"
        missing_deps=1
    fi
    
    if [[ $missing_deps -eq 1 ]]; then
        echo -e "${RED}Please ensure all required scripts are present and executable${NC}"
        exit 1
    fi
}

# Display configuration and job discovery
show_info() {
    if [[ $QUIET -eq 0 ]]; then
        echo -e "${BLUE}AlphaFold3 Output Sync Orchestrator${NC}"
        echo "===================================="
        echo "Base directory: $BASE_DIR"
        echo "Configuration: $CONFIG_FILE"
        echo
        
        if [[ $DRY_RUN -eq 1 ]]; then
            echo -e "${YELLOW}DRY RUN MODE - No jobs will be submitted${NC}"
            echo
        fi
        
        # Show job discovery info
        "$JOB_DISCOVERY" info
        echo
    fi
}

# Create output directories
create_directories() {
    if [[ $DRY_RUN -eq 0 ]]; then
        mkdir -p "${OUTPUT_DIR}"
        mkdir -p "${OUTPUT_DIR}/msa"
        mkdir -p "${LOGS_DIR}"
    fi
}

# Submit rsync job
submit_rsync_job() {
    if [[ $SKIP_RSYNC -eq 1 ]]; then
        echo -e "${YELLOW}Skipping rsync job submission${NC}"
        return 0
    fi
    
    echo -e "${GREEN}Step 1: Preparing parallel rsync job${NC}"
    
    # Get job chunks dynamically
    local chunk_info
    if ! chunk_info=$("$JOB_DISCOVERY" chunks 2>/dev/null); then
        echo -e "${RED}Error: Failed to discover jobs for rsync${NC}"
        return 1
    fi
    
    # Parse chunk information
    eval "$chunk_info"
    
    if [[ $TOTAL_JOBS -eq 0 ]]; then
        echo -e "${YELLOW}No jobs found to sync${NC}"
        return 0
    fi
    
    echo "Jobs discovered: $TOTAL_JOBS"
    echo "Array size: $ARRAY_SIZE"
    echo "Jobs per task: $JOBS_PER_TASK"
    echo "Total parallelism: $(( ARRAY_SIZE * PARALLEL_WORKERS )) workers"
    
    if [[ $DRY_RUN -eq 1 ]]; then
        echo -e "${YELLOW}Would submit: sbatch --array=1-${ARRAY_SIZE}%${MAX_ARRAY_TASKS} $RSYNC_SBATCH${NC}"
        return 0
    fi
    
    # Export environment variables for the array job
    export SYNC_JOBS_LIST="$JOBS_LIST"
    export SYNC_TOTAL_JOBS="$TOTAL_JOBS"
    export SYNC_JOBS_PER_TASK="$JOBS_PER_TASK"
    
    # Submit the rsync array job
    local rsync_job_id
    if rsync_job_id=$(sbatch --parsable --array=1-${ARRAY_SIZE}%${MAX_ARRAY_TASKS} \
                            --cpus-per-task="$CPUS_PER_TASK" \
                            --mem="$MEMORY_PER_TASK" \
                            --time="$TIME_LIMIT" \
                            --partition="$PARTITION" \
                            "$RSYNC_SBATCH"); then
        echo "Rsync job submitted: ID $rsync_job_id"
        echo "Monitor with: squeue -j $rsync_job_id"
        echo
        
        # Store job ID for dependency management
        RSYNC_JOB_ID="$rsync_job_id"
        return 0
    else
        echo -e "${RED}Failed to submit rsync job${NC}"
        return 1
    fi
}

# Count total job groups for MSA archiving
count_msa_job_groups() {
    # Fast job group counting using the same logic as the sbatch script
    local temp_file=$(mktemp)
    
    # Main jobs directory
    find "$JOBS_DIR" -maxdepth 1 -type d -name "*-*" -exec basename {} \; | grep -v "human_test_set" | cut -d'-' -f1 >> "$temp_file"
    
    # Human test set subdirectory
    if [[ -d "$JOBS_DIR/human_test_set" ]]; then
        find "$JOBS_DIR/human_test_set" -maxdepth 1 -type d -name "*-*" -exec basename {} \; | cut -d'-' -f1 >> "$temp_file"
    fi
    
    # Count unique groups
    local count=$(sort -u "$temp_file" | wc -l)
    rm -f "$temp_file"
    echo $count
}

# Submit MSA archiving job
submit_msa_job() {
    if [[ $SKIP_MSA -eq 1 ]]; then
        echo -e "${YELLOW}Skipping MSA archiving job submission${NC}"
        return 0
    fi
    
    echo -e "${GREEN}Step 2: Submitting MSA archiving job${NC}"
    
    # Calculate total number of job groups for dynamic array sizing
    local total_msa_groups
    echo "Discovering job groups for MSA archiving..."
    total_msa_groups=$(count_msa_job_groups)
    echo "Found $total_msa_groups job groups to process"
    
    if [[ $total_msa_groups -eq 0 ]]; then
        echo -e "${YELLOW}No job groups found - skipping MSA archiving${NC}"
        return 0
    fi
    
    # Use same concurrency limit as rsync (from config)
    local max_concurrent=${MAX_ARRAY_TASKS:-10}
    
    if [[ $DRY_RUN -eq 1 ]]; then
        echo -e "${YELLOW}Would submit: sbatch --array=1-${total_msa_groups}%${max_concurrent} $MSA_SBATCH${NC}"
        return 0
    fi
    
    # Submit MSA archiving with dynamic array sizing
    local msa_job_id
    if msa_job_id=$(sbatch --parsable --array=1-${total_msa_groups}%${max_concurrent} \
                          --cpus-per-task=4 \
                          --mem=8G \
                          --time=6:00:00 \
                          --partition="${PARTITION}" \
                          "$MSA_SBATCH"); then
        echo "MSA archiving job submitted: ID $msa_job_id"
        echo "Array size: 1-${total_msa_groups}%${max_concurrent}"
        echo "Monitor with: squeue -j $msa_job_id"
        echo
        return 0
    else
        echo -e "${RED}Failed to submit MSA archiving job${NC}"
        return 1
    fi
}

# Submit seed compression job
submit_compress_job() {
    if [[ $SKIP_COMPRESS -eq 1 ]]; then
        echo -e "${YELLOW}Skipping seed compression job submission${NC}"
        return 0
    fi
    
    echo -e "${GREEN}Step 3: Submitting seed compression job${NC}"
    
    local dependency_flag=""
    if [[ -n "${RSYNC_JOB_ID:-}" ]] && [[ $DRY_RUN -eq 0 ]]; then
        dependency_flag="--dependency=afterany:$RSYNC_JOB_ID"
        echo "Setting dependency on rsync job: $RSYNC_JOB_ID"
    fi
    
    if [[ $DRY_RUN -eq 1 ]]; then
        echo -e "${YELLOW}Would submit: sbatch $dependency_flag $COMPRESS_SBATCH${NC}"
        return 0
    fi
    
    # Submit seed compression (depends on rsync completion)
    local compress_job_id
    if compress_job_id=$(sbatch --parsable $dependency_flag "$COMPRESS_SBATCH"); then
        echo "Seed compression job submitted: ID $compress_job_id"
        echo "Monitor with: squeue -j $compress_job_id"
        echo
        return 0
    else
        echo -e "${RED}Failed to submit seed compression job${NC}"
        return 1
    fi
}

# Show monitoring commands
show_monitoring_info() {
    if [[ $DRY_RUN -eq 0 && $QUIET -eq 0 ]]; then
        echo -e "${GREEN}Job Monitoring${NC}"
        echo "=============="
        echo "Check all jobs: squeue -u \$USER"
        echo "View logs: ls -la logs/sync_rsync_*"
        echo "           ls -la logs/archive_msa_*"
        echo "           ls -la logs/pack-seeds_*"
        echo
        echo "Progress tracking:"
        echo "  tail -f logs/sync_rsync_*.out     # Rsync progress"
        echo "  tail -f logs/archive_msa_*.out    # MSA archiving"
        echo "  tail -f logs/pack-seeds_*.out     # Seed compression"
    fi
}

# Main execution
main() {
    check_dependencies
    show_info
    create_directories
    
    local exit_code=0
    
    # Submit jobs in order
    if ! submit_rsync_job; then
        exit_code=1
    fi
    
    if ! submit_msa_job; then
        exit_code=1
    fi
    
    if ! submit_compress_job; then
        exit_code=1
    fi
    
    # Summary
    if [[ $DRY_RUN -eq 0 ]]; then
        if [[ $exit_code -eq 0 ]]; then
            echo -e "${GREEN}All jobs submitted successfully!${NC}"
            echo
            show_monitoring_info
        else
            echo -e "${RED}Some jobs failed to submit${NC}"
            return $exit_code
        fi
    else
        echo -e "${YELLOW}Dry run complete - no jobs submitted${NC}"
        echo "Run without -n flag to submit jobs"
    fi
    
    return $exit_code
}

# Execute main function
main "$@"