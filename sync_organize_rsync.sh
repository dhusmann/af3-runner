#!/bin/bash

# sync_organize_rsync.sh - Parallel rsync worker for AlphaFold3 output organization
#
# This script performs the actual rsync operations with parallel processing support.
# It maintains all the organization logic from the original sync_organize_outputs.sh
# but processes jobs in parallel using GNU parallel.
#
# Usage:
#   echo "job_list" | ./sync_organize_rsync.sh [mode] [parallel_workers]
#   ./sync_organize_rsync.sh process_file job_list.txt [parallel_workers]
#
# Modes:
#   parallel   - Use GNU parallel (default)
#   sequential - Process one at a time
#   process_file - Read from file instead of stdin

set -euo pipefail

# Load configuration
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
CONFIG_FILE="${SCRIPT_DIR}/sync_parallel.conf"

if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    # Default settings if config not found
    BASE_DIR="/scratch/groups/ogozani/alphafold3"
    OUTPUT_DIR="${BASE_DIR}/output"
    PARALLEL_WORKERS=16
    RSYNC_TIMEOUT=60
    PROGRESS_REPORT_INTERVAL=100
    DEBUG_MODE=false
    RSYNC_EXTRA_OPTS="--sparse --inplace"
fi

# Command line arguments
MODE="${1:-parallel}"
# Default: assume workers argument immediately follows mode unless process_file
WORKERS="${2:-$PARALLEL_WORKERS}"

# For process_file mode we expect: process_file <jobs_file> [workers]
if [[ "$MODE" == "process_file" ]]; then
    WORKERS="${3:-$PARALLEL_WORKERS}"
fi

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Statistics
SYNC_COUNT=0
SKIP_COUNT=0
ERROR_COUNT=0

# Ensure output directories exist
create_output_directories() {
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}/msa"
    mkdir -p "${OUTPUT_DIR}/human_test_set"
    mkdir -p "${OUTPUT_DIR}/arabidopsis_EEF1A"
    mkdir -p "${OUTPUT_DIR}/ecoli/ecEFTU"
    mkdir -p "${OUTPUT_DIR}/ecoli/ecRPL11"
    mkdir -p "${OUTPUT_DIR}/dpEEF1A"
    mkdir -p "${OUTPUT_DIR}/RPL29"
    mkdir -p "${OUTPUT_DIR}/RPL36A"
    mkdir -p "${OUTPUT_DIR}/SETD6"
    mkdir -p "${OUTPUT_DIR}/human_EF_KMTs/EEF1A"
    mkdir -p "${OUTPUT_DIR}/human_EF_KMTs/EF2"
    mkdir -p "${OUTPUT_DIR}/yeast/EFM"
    mkdir -p "${OUTPUT_DIR}/yeast/RKM"
    mkdir -p "${OUTPUT_DIR}/yeast/SET"
    mkdir -p "${OUTPUT_DIR}/NSD2i"
}

# Function to determine destination based on job name (from original script)
get_destination() {
    local job_name="$1"
    local job_lower=$(echo "$job_name" | tr '[:upper:]' '[:lower:]')
    
    # Check patterns in order of preference (same as original)
    if [[ "$job_lower" == *"ateef1a"* ]]; then
        echo "arabidopsis_EEF1A"
    elif [[ "$job_lower" == *"eftu"* ]]; then
        echo "ecoli/ecEFTU"
    elif [[ "$job_lower" == *"ecrpl"* ]]; then
        echo "ecoli/ecRPL11"
    elif [[ "$job_lower" == *"dpeef1a"* ]]; then
        echo "dpEEF1A"
    elif [[ "$job_lower" == *"rpl29"* ]]; then
        echo "RPL29"
    elif [[ "$job_lower" == *"rpl36a"* ]]; then
        echo "RPL36A"
    elif [[ "$job_lower" == *"setd6"* ]]; then
        echo "SETD6"
    elif [[ "$job_lower" =~ ^eef1akmt.* ]] || [[ "$job_lower" =~ ^mettl13.* ]]; then
        echo "human_EF_KMTs/EEF1A"
    elif [[ "$job_lower" =~ ^fam86.*eef.* ]]; then
        echo "human_EF_KMTs/EF2"
    elif [[ "$job_lower" =~ ^scefm.* ]]; then
        echo "yeast/EFM"
    elif [[ "$job_lower" =~ ^scrkm.* ]]; then
        echo "yeast/RKM"
    elif [[ "$job_lower" =~ ^spset.* ]]; then
        echo "yeast/SET"
    elif [[ "$job_lower" =~ .*aa.*[0-9]+to[0-9]+.* ]]; then
        echo "NSD2i"
    else
        echo ""  # No match - will go to top level
    fi
}

# Function to process a single job
process_single_job() {
    local job_line="$1"
    local job_type=$(echo "$job_line" | cut -d'|' -f1)
    local job_dir=$(echo "$job_line" | cut -d'|' -f2)
    local output_subdir=$(echo "$job_line" | cut -d'|' -f3)
    
    local job_name=$(basename "$job_dir")
    local subdir_name=$(basename "$output_subdir")
    
    if [[ $DEBUG_MODE == "true" ]]; then
        echo "Processing: $job_type | $job_name | $subdir_name" >&2
    fi
    
    # Skip if source doesn't exist
    if [[ ! -d "$output_subdir" ]]; then
        echo "SKIP|$job_name|Source not found: $output_subdir" >&2
        return 0
    fi
    
    local dest_path=""
    local final_destination=""
    
    if [[ "$job_type" == "HUMAN_TEST_SET" ]]; then
        # Human test set jobs always go to human_test_set/
        final_destination="${OUTPUT_DIR}/human_test_set/${subdir_name}/"
    else
        # Regular jobs use pattern matching
        dest_path=$(get_destination "$subdir_name")
        
        if [[ -n "$dest_path" ]]; then
            final_destination="${OUTPUT_DIR}/${dest_path}/${subdir_name}/"
        else
            # No pattern match - goes to top level
            final_destination="${OUTPUT_DIR}/${subdir_name}/"
        fi
    fi
    
    # Ensure destination directory exists
    mkdir -p "$(dirname "$final_destination")"
    
    # Build rsync command
    local rsync_cmd=(
        rsync
        -a
        --exclude="*_data.json"
        --exclude="TERMS_OF_USE.md"
        --timeout="$RSYNC_TIMEOUT"
    )
    
    # Add extra options if configured
    if [[ -n "$RSYNC_EXTRA_OPTS" ]]; then
        read -ra extra_opts <<< "$RSYNC_EXTRA_OPTS"
        rsync_cmd+=("${extra_opts[@]}")
    fi
    
    rsync_cmd+=("$output_subdir/" "$final_destination")
    
    # Execute rsync
    if "${rsync_cmd[@]}" 2>/dev/null; then
        echo "SUCCESS|$job_name|$subdir_name|$dest_path"
        return 0
    else
        echo "ERROR|$job_name|$subdir_name|rsync failed" >&2
        return 1
    fi
}

# Export the function for GNU parallel
export -f process_single_job
export -f get_destination
export OUTPUT_DIR RSYNC_TIMEOUT RSYNC_EXTRA_OPTS DEBUG_MODE

# Function to process jobs in parallel
process_parallel() {
    local input_source="$1"
    
    echo "Starting parallel processing with $WORKERS workers..." >&2
    
    if [[ "$input_source" == "stdin" ]]; then
        parallel -j "$WORKERS" --line-buffer process_single_job
    else
        parallel -j "$WORKERS" --line-buffer -a "$input_source" process_single_job
    fi
}

# Function to process jobs sequentially
process_sequential() {
    local input_source="$1"
    local line_count=0
    
    echo "Starting sequential processing..." >&2
    
    if [[ "$input_source" == "stdin" ]]; then
        while IFS= read -r job_line; do
            ((line_count++))
            if [[ $((line_count % PROGRESS_REPORT_INTERVAL)) -eq 0 ]]; then
                echo "Processed $line_count jobs..." >&2
            fi
            process_single_job "$job_line"
        done
    else
        while IFS= read -r job_line; do
            ((line_count++))
            if [[ $((line_count % PROGRESS_REPORT_INTERVAL)) -eq 0 ]]; then
                echo "Processed $line_count jobs..." >&2
            fi
            process_single_job "$job_line"
        done < "$input_source"
    fi
}

# Function to collect and report statistics
collect_statistics() {
    local tmp_file
    tmp_file=$(mktemp)
    cat > "$tmp_file"

    local success_count=$(grep -c "^SUCCESS" "$tmp_file" || true)
    local skip_count=$(grep -c "^SKIP" "$tmp_file" || true) 
    local error_count=$(grep -c "^ERROR" "$tmp_file" || true)
    local total_count=$((success_count + skip_count + error_count))
    
    echo >&2
    echo "Processing Statistics:" >&2
    echo "  Total processed: $total_count" >&2
    echo "  Successful syncs: $success_count" >&2
    echo "  Skipped (missing): $skip_count" >&2
    echo "  Errors: $error_count" >&2
    
    if [[ $error_count -gt 0 ]]; then
        echo -e "${RED}Warning: $error_count jobs had errors${NC}" >&2
        echo "Check error messages above for details" >&2
        rm -f "$tmp_file"
        return 1
    fi
    
    rm -f "$tmp_file"
    return 0
}

# Main execution
main() {
    echo -e "${BLUE}AlphaFold3 Parallel Rsync Worker${NC}" >&2
    echo "=================================" >&2
    echo "Mode: $MODE" >&2
    echo "Workers: $WORKERS" >&2
    echo "Output directory: $OUTPUT_DIR" >&2
    echo >&2
    
    # Create output directories
    create_output_directories
    
    # Determine input source
    local input_source="stdin"
    if [[ "$MODE" == "process_file" ]]; then
        if [[ $# -lt 2 ]]; then
            echo "Error: process_file mode requires filename argument" >&2
            exit 1
        fi
        input_source="$2"
        if [[ ! -f "$input_source" ]]; then
            echo "Error: Input file not found: $input_source" >&2
            exit 1
        fi
    fi
    
    # Check if GNU parallel is available for parallel mode
    if [[ "$MODE" == "parallel" ]] && ! command -v parallel >/dev/null 2>&1; then
        echo "Warning: GNU parallel not found, falling back to sequential mode" >&2
        MODE="sequential"
    fi
    
    # Process based on mode
    case "$MODE" in
        "parallel")
            process_parallel "$input_source" | collect_statistics
            ;;
        "sequential")
            process_sequential "$input_source" | collect_statistics
            ;;
        "process_file")
            process_parallel "$input_source" | collect_statistics
            ;;
        *)
            echo "Unknown mode: $MODE" >&2
            echo "Valid modes: parallel, sequential, process_file" >&2
            exit 1
            ;;
    esac
}

# Check if being run directly (not sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
