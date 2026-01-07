#!/bin/bash

# get_job_status_detailed.sh - Extended job status report with export options
# Additional features: export to CSV, filter by stage, check specific seeds

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Default values
BASE_DIR="/scratch/groups/ogozani/alphafold3"
CSV_FILE=""
EXPORT_FILE=""
FILTER_STAGE=""
CHECK_SEED=""
SHOW_PATHS=0

# Usage function
usage() {
    echo "Usage: $0 [-f CSV_FILE] [-e EXPORT_FILE] [-stage N] [-seed N] [-p] [-h]"
    echo "  -f CSV_FILE      Check only jobs listed in CSV file"
    echo "  -e EXPORT_FILE   Export results to CSV file"
    echo "  -stage N         Filter by stage (1, 2, or 3)"
    echo "  -seed N          Check which jobs have completed seed N"
    echo "  -p               Show full paths to job directories"
    echo "  -h               Show this help message"
    echo
    echo "Examples:"
    echo "  $0 -stage 2                    # Show only Stage 2 jobs"
    echo "  $0 -e status_report.csv         # Export all job status to CSV"
    echo "  $0 -seed 0 -stage 3             # Show Stage 3 jobs with seed 0 complete"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -f)
            CSV_FILE="$2"
            if [ ! -f "$CSV_FILE" ]; then
                echo "Error: CSV file not found: $CSV_FILE"
                exit 1
            fi
            shift 2
            ;;
        -e)
            EXPORT_FILE="$2"
            shift 2
            ;;
        -stage)
            FILTER_STAGE="$2"
            if [[ ! "$FILTER_STAGE" =~ ^[123]$ ]]; then
                echo "Error: Stage must be 1, 2, or 3"
                exit 1
            fi
            shift 2
            ;;
        -seed)
            CHECK_SEED="$2"
            if [[ ! "$CHECK_SEED" =~ ^[0-9]+$ ]]; then
                echo "Error: Seed must be a number"
                exit 1
            fi
            shift 2
            ;;
        -p)
            SHOW_PATHS=1
            shift
            ;;
        -h)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Function to find job directory
find_job_directory() {
    local job_name=$1
    
    if [ -d "${BASE_DIR}/jobs/$job_name" ]; then
        echo "${BASE_DIR}/jobs/$job_name"
        return 0
    elif [ -d "${BASE_DIR}/jobs/human_test_set/$job_name" ]; then
        echo "${BASE_DIR}/jobs/human_test_set/$job_name"
        return 0
    fi
    
    return 1
}

# Function to check job stage with detailed info
check_job_stage_detailed() {
    local job_name=$1
    local job_dir=$2
    local result=""
    
    # Check for alphafold_input.json
    if [ ! -f "${job_dir}/alphafold_input.json" ]; then
        echo "0|No input JSON"
        return
    fi
    
    # Check for MSA output
    local msa_json=$(find "${job_dir}/output_msa" -name "*.json" -type f 2>/dev/null | head -1)
    if [ -z "$msa_json" ]; then
        echo "1|Ready for MSA"
        return
    fi
    
    # Check for GPU output
    if [ ! -d "${job_dir}/output" ] || [ -z "$(ls -A "${job_dir}/output" 2>/dev/null)" ]; then
        local msa_date=$(stat -c "%Y" "$msa_json" 2>/dev/null)
        local msa_human=$(date -d "@$msa_date" "+%Y-%m-%d %H:%M" 2>/dev/null || echo "Unknown")
        echo "2|MSA complete at $msa_human"
        return
    fi
    
    # Get completion time from newest file in output
    local newest=$(find "${job_dir}/output" -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2-)
    local gpu_date=""
    if [ -n "$newest" ]; then
        gpu_date=$(stat -c "%Y" "$newest" 2>/dev/null)
        gpu_human=$(date -d "@$gpu_date" "+%Y-%m-%d %H:%M" 2>/dev/null || echo "Unknown")
    else
        gpu_human="Unknown"
    fi
    
    echo "3|GPU complete at $gpu_human"
}

# Function to get completed seeds with details
get_completed_seeds_detailed() {
    local job_dir=$1
    local check_specific=$2
    local seeds=()
    local has_specific=0
    
    # Look for subdirectories in output/
    for subdir in "${job_dir}/output/"*/; do
        if [ -d "$subdir" ]; then
            # Look for seed directories within the subdirectory
            for seed_dir in "${subdir}"seed-*_sample-*/; do
                if [ -d "$seed_dir" ]; then
                    # Extract seed number from directory name
                    seed_num=$(basename "$seed_dir" | grep -oP 'seed-\K\d+(?=_sample)')
                    
                    # Check if this seed is already in our list
                    if [[ ! " ${seeds[@]} " =~ " ${seed_num} " ]]; then
                        # Check if all 5 samples (0-4) have model.cif for this seed
                        local all_samples_complete=1
                        for sample in 0 1 2 3 4; do
                            if [ ! -f "${subdir}seed-${seed_num}_sample-${sample}/model.cif" ]; then
                                all_samples_complete=0
                                break
                            fi
                        done
                        
                        if [ $all_samples_complete -eq 1 ]; then
                            seeds+=($seed_num)
                            if [ -n "$check_specific" ] && [ "$seed_num" -eq "$check_specific" ]; then
                                has_specific=1
                            fi
                        fi
                    fi
                fi
            done
        fi
    done
    
    # Return results
    if [ -n "$check_specific" ]; then
        echo "$has_specific"
    else
        if [ ${#seeds[@]} -gt 0 ]; then
            printf '%s\n' "${seeds[@]}" | sort -n | tr '\n' ',' | sed 's/,$//'
        fi
    fi
}

# Initialize export file if requested
if [ -n "$EXPORT_FILE" ]; then
    echo "job_name,stage,status,completed_seeds,path" > "$EXPORT_FILE"
fi

# Header
echo "AlphaFold3 Job Status Report"
echo "============================"
if [ -n "$FILTER_STAGE" ]; then
    echo -e "${CYAN}Filter: Stage $FILTER_STAGE${NC}"
fi
if [ -n "$CHECK_SEED" ]; then
    echo -e "${CYAN}Filter: Jobs with seed $CHECK_SEED complete${NC}"
fi
echo

# Get list of jobs to process
if [ -n "$CSV_FILE" ]; then
    echo "Reading jobs from: $CSV_FILE"
    JOBS=$(tail -n +2 "$CSV_FILE" | grep -v "^$" | tr -d '\r' | awk '{print $1}')
else
    echo "Scanning all job directories..."
    JOBS=""
    
    # Scan main jobs directory
    if [ -d "${BASE_DIR}/jobs" ]; then
        for job_dir in "${BASE_DIR}/jobs/"*/; do
            if [ -d "$job_dir" ]; then
                job_name=$(basename "$job_dir")
                if [[ ! "$job_name" =~ ^human_test_set$ ]]; then
                    JOBS="$JOBS$job_name"$'\n'
                fi
            fi
        done
    fi
    
    # Scan human_test_set
    if [ -d "${BASE_DIR}/jobs/human_test_set" ]; then
        for job_dir in "${BASE_DIR}/jobs/human_test_set/"*/; do
            if [ -d "$job_dir" ]; then
                job_name=$(basename "$job_dir")
                JOBS="$JOBS$job_name"$'\n'
            fi
        done
    fi
    
    JOBS=$(echo "$JOBS" | grep -v "^$" | sort)
fi

# Process each job
TOTAL_JOBS=0
DISPLAYED_JOBS=0
STAGE_COUNTS=(0 0 0 0)  # Not found, Stage 1, Stage 2, Stage 3

echo "Job Status:"
echo "-----------"

while IFS= read -r job_name; do
    if [ -n "$job_name" ]; then
        ((TOTAL_JOBS++))
        
        job_dir=$(find_job_directory "$job_name")
        
        if [ -z "$job_dir" ]; then
            if [ -z "$FILTER_STAGE" ]; then
                echo -e "${RED}[NOT FOUND]${NC} $job_name"
                ((DISPLAYED_JOBS++))
            fi
            ((STAGE_COUNTS[0]++))
            continue
        fi
        
        # Get detailed stage info
        stage_info=$(check_job_stage_detailed "$job_name" "$job_dir")
        stage=$(echo "$stage_info" | cut -d'|' -f1)
        status=$(echo "$stage_info" | cut -d'|' -f2)
        
        # Apply stage filter
        if [ -n "$FILTER_STAGE" ] && [ "$stage" != "$FILTER_STAGE" ]; then
            ((STAGE_COUNTS[$stage]++))
            continue
        fi
        
        # Get seed information for stage 3
        seeds=""
        if [ "$stage" -eq 3 ]; then
            if [ -n "$CHECK_SEED" ]; then
                has_seed=$(get_completed_seeds_detailed "$job_dir" "$CHECK_SEED")
                if [ "$has_seed" -eq 0 ]; then
                    ((STAGE_COUNTS[3]++))
                    continue
                fi
            fi
            seeds=$(get_completed_seeds_detailed "$job_dir")
        fi
        
        # Format output
        case $stage in
            1)
                output="${YELLOW}[STAGE 1]${NC} $job_name - $status"
                ((STAGE_COUNTS[1]++))
                ;;
            2)
                output="${BLUE}[STAGE 2]${NC} $job_name - $status"
                ((STAGE_COUNTS[2]++))
                ;;
            3)
                if [ -n "$seeds" ]; then
                    output="${GREEN}[STAGE 3]${NC} $job_name - $status | Seeds: $seeds"
                else
                    output="${GREEN}[STAGE 3]${NC} $job_name - $status | No complete seeds"
                fi
                ((STAGE_COUNTS[3]++))
                ;;
        esac
        
        # Add path if requested
        if [ $SHOW_PATHS -eq 1 ]; then
            output="$output\n         Path: $job_dir"
        fi
        
        echo -e "$output"
        ((DISPLAYED_JOBS++))
        
        # Export to CSV if requested
        if [ -n "$EXPORT_FILE" ]; then
            echo "$job_name,$stage,$status,\"$seeds\",\"$job_dir\"" >> "$EXPORT_FILE"
        fi
    fi
done <<< "$JOBS"

# Summary
echo
echo "Summary"
echo "======="
echo "Total jobs checked:   $TOTAL_JOBS"
echo "Jobs displayed:       $DISPLAYED_JOBS"
echo
echo -e "${YELLOW}Stage 1 (Need MSA):   ${STAGE_COUNTS[1]}${NC}"
echo -e "${BLUE}Stage 2 (Need GPU):   ${STAGE_COUNTS[2]}${NC}"
echo -e "${GREEN}Stage 3 (Complete):   ${STAGE_COUNTS[3]}${NC}"
if [ ${STAGE_COUNTS[0]} -gt 0 ]; then
    echo -e "${RED}Not found:            ${STAGE_COUNTS[0]}${NC}"
fi

# Export notification
if [ -n "$EXPORT_FILE" ]; then
    echo
    echo -e "${CYAN}Results exported to: $EXPORT_FILE${NC}"
fi

echo
echo "Done."
