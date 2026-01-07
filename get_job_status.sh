#!/bin/bash

# get_job_status.sh - Report status of AlphaFold3 jobs
# Usage: 
#   ./get_job_status.sh                    # Check all jobs in jobs/ directories
#   ./get_job_status.sh -f jobs.csv        # Check specific jobs from CSV
#   ./get_job_status.sh -s                # Summary only
#   ./get_job_status.sh -v                # Verbose output

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
BASE_DIR="/scratch/groups/ogozani/alphafold3"
CSV_FILE=""
SUMMARY_ONLY=0
VERBOSE=0

# Usage function
usage() {
    echo "Usage: $0 [-f CSV_FILE] [-s] [-v] [-h]"
    echo "  -f CSV_FILE    Check only jobs listed in CSV file"
    echo "  -s             Summary only (no individual job listing)"
    echo "  -v             Verbose output"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "f:svh" opt; do
    case $opt in
        f)
            CSV_FILE=$OPTARG
            if [ ! -f "$CSV_FILE" ]; then
                echo "Error: CSV file not found: $CSV_FILE"
                exit 1
            fi
            ;;
        s)
            SUMMARY_ONLY=1
            ;;
        v)
            VERBOSE=1
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
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

# Function to check job stage
check_job_stage() {
    local job_name=$1
    local job_dir=$2
    
    # Check for alphafold_input.json
    if [ ! -f "${job_dir}/alphafold_input.json" ]; then
        echo "0"  # No input file
        return
    fi
    
    # Check for MSA output
    local msa_json=$(find "${job_dir}/output_msa" -name "*.json" -type f 2>/dev/null | head -1)
    if [ -z "$msa_json" ]; then
        echo "1"  # Stage 1: Ready for MSA
        return
    fi
    
    # Check for GPU output
    if [ ! -d "${job_dir}/output" ] || [ -z "$(ls -A "${job_dir}/output" 2>/dev/null)" ]; then
        echo "2"  # Stage 2: MSA complete, ready for GPU
        return
    fi
    
    echo "3"  # Stage 3: GPU complete
}

# Function to get completed seeds for stage 3 jobs
get_completed_seeds() {
    local job_dir=$1
    local seeds=()
    
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
                        fi
                    fi
                fi
            done
        fi
    done
    
    # Sort seeds numerically
    if [ ${#seeds[@]} -gt 0 ]; then
        printf '%s\n' "${seeds[@]}" | sort -n | tr '\n' ',' | sed 's/,$//'
    fi
}

# Function to process a single job
process_job() {
    local job_name=$1
    local job_dir=$(find_job_directory "$job_name")
    
    if [ -z "$job_dir" ]; then
        if [ $VERBOSE -eq 1 ]; then
            echo -e "${RED}[NOT FOUND]${NC} $job_name"
        fi
        return 1
    fi
    
    local stage=$(check_job_stage "$job_name" "$job_dir")
    
    case $stage in
        0)
            if [ $VERBOSE -eq 1 ]; then
                echo -e "${RED}[NO INPUT]${NC} $job_name"
            fi
            return 2
            ;;
        1)
            echo -e "${YELLOW}[STAGE 1]${NC} $job_name"
            return 3
            ;;
        2)
            echo -e "${BLUE}[STAGE 2]${NC} $job_name"
            return 4
            ;;
        3)
            local seeds=$(get_completed_seeds "$job_dir")
            if [ -n "$seeds" ]; then
                echo -e "${GREEN}[STAGE 3]${NC} $job_name | Seeds: $seeds"
            else
                echo -e "${GREEN}[STAGE 3]${NC} $job_name | No complete seeds found"
            fi
            return 5
            ;;
    esac
}

# Initialize counters
TOTAL_JOBS=0
NOT_FOUND=0
NO_INPUT=0
STAGE_1=0
STAGE_2=0
STAGE_3=0

# Main processing
echo "AlphaFold3 Job Status Report"
echo "============================"
echo

# Get list of jobs to process
if [ -n "$CSV_FILE" ]; then
    # Read jobs from CSV
    echo "Reading jobs from: $CSV_FILE"
    
    # Skip header and empty lines
    JOBS=$(tail -n +2 "$CSV_FILE" | grep -v "^$" | tr -d '\r' | awk '{print $1}')
else
    # Find all jobs in both directories
    echo "Scanning all job directories..."
    
    JOBS=""
    
    # Scan main jobs directory
    if [ -d "${BASE_DIR}/jobs" ]; then
        for job_dir in "${BASE_DIR}/jobs/"*/; do
            if [ -d "$job_dir" ]; then
                job_name=$(basename "$job_dir")
                # Skip subdirectories that might be category folders
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
    
    # Remove empty lines and sort
    JOBS=$(echo "$JOBS" | grep -v "^$" | sort)
fi

# Count total jobs
TOTAL_JOBS=$(echo "$JOBS" | grep -v "^$" | wc -l)
echo "Total jobs to check: $TOTAL_JOBS"
echo

# Process each job
if [ $SUMMARY_ONLY -eq 0 ]; then
    echo "Individual Job Status:"
    echo "---------------------"
fi

while IFS= read -r job_name; do
    if [ -n "$job_name" ]; then
        if [ $SUMMARY_ONLY -eq 0 ]; then
            process_job "$job_name"
        else
            # Just count for summary
            job_dir=$(find_job_directory "$job_name")
            if [ -z "$job_dir" ]; then
                ((NOT_FOUND++))
                continue
            fi
            
            stage=$(check_job_stage "$job_name" "$job_dir")
            case $stage in
                0) ((NO_INPUT++)) ;;
                1) ((STAGE_1++)) ;;
                2) ((STAGE_2++)) ;;
                3) ((STAGE_3++)) ;;
            esac
        fi
        
        # Update counters based on return value
        case $? in
            1) ((NOT_FOUND++)) ;;
            2) ((NO_INPUT++)) ;;
            3) ((STAGE_1++)) ;;
            4) ((STAGE_2++)) ;;
            5) ((STAGE_3++)) ;;
        esac
    fi
done <<< "$JOBS"

# Summary
echo
echo "Summary"
echo "======="
echo "Total jobs:           $TOTAL_JOBS"
if [ $NOT_FOUND -gt 0 ] || [ $NO_INPUT -gt 0 ]; then
    echo -e "${RED}Not found:            $NOT_FOUND${NC}"
    echo -e "${RED}No input JSON:        $NO_INPUT${NC}"
fi
echo -e "${YELLOW}Stage 1 (Need MSA):   $STAGE_1${NC}"
echo -e "${BLUE}Stage 2 (Need GPU):   $STAGE_2${NC}"
echo -e "${GREEN}Stage 3 (Complete):   $STAGE_3${NC}"

# Calculate percentage complete
if [ $TOTAL_JOBS -gt 0 ]; then
    PERCENT_COMPLETE=$(( (STAGE_3 * 100) / TOTAL_JOBS ))
    echo
    echo "Progress: ${PERCENT_COMPLETE}% complete"
fi

# Additional statistics if verbose
if [ $VERBOSE -eq 1 ] && [ $STAGE_3 -gt 0 ]; then
    echo
    echo "Checking seed completion for Stage 3 jobs..."
    TOTAL_SEEDS=0
    JOBS_WITH_SEEDS=0
    
    while IFS= read -r job_name; do
        if [ -n "$job_name" ]; then
            job_dir=$(find_job_directory "$job_name")
            if [ -n "$job_dir" ]; then
                stage=$(check_job_stage "$job_name" "$job_dir")
                if [ "$stage" -eq 3 ]; then
                    seeds=$(get_completed_seeds "$job_dir")
                    if [ -n "$seeds" ]; then
                        seed_count=$(echo "$seeds" | tr ',' '\n' | wc -l)
                        TOTAL_SEEDS=$((TOTAL_SEEDS + seed_count))
                        ((JOBS_WITH_SEEDS++))
                        
                        if [ $VERBOSE -eq 1 ]; then
                            echo "  $job_name: $seed_count seeds ($seeds)"
                        fi
                    fi
                fi
            fi
        fi
    done <<< "$JOBS"
    
    if [ $JOBS_WITH_SEEDS -gt 0 ]; then
        AVG_SEEDS=$(( TOTAL_SEEDS / JOBS_WITH_SEEDS ))
        echo
        echo "Stage 3 jobs with completed seeds: $JOBS_WITH_SEEDS / $STAGE_3"
        echo "Total completed seeds: $TOTAL_SEEDS"
        echo "Average seeds per job: $AVG_SEEDS"
    else
        echo "No Stage 3 jobs have completed seeds (check output directory structure)"
    fi
fi

# Check for running jobs
echo
echo "Currently Running Jobs:"
echo "----------------------"
MSA_RUNNING=$(squeue -u $USER -n af3_msa_array -h | wc -l)
GPU_RUNNING=$(squeue -u $USER -p gpu -h | wc -l)

echo "MSA array jobs running: $MSA_RUNNING"
echo "GPU jobs running:       $GPU_RUNNING"

echo
echo "Done."
