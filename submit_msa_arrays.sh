#!/bin/bash

# Script to submit MSA jobs as SLURM arrays from CSV file
# Splits jobs between hns and normal partitions for better parallelization

# Configuration
AF3_BASE_DIR="/scratch/groups/ogozani/alphafold3"
CSV_FILE="${1:-msa_array_jobs.csv}"
MAX_ARRAY_SIZE=1000
SCRIPT_DIR="$(dirname "$0")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}AlphaFold 3 MSA Array Job Submission (Multi-Partition)${NC}"
echo "====================================================="

# Check if CSV file exists
if [ ! -f "$CSV_FILE" ]; then
    echo -e "${RED}Error: CSV file not found: $CSV_FILE${NC}"
    exit 1
fi

# Check if array submission script exists
ARRAY_SCRIPT="${SCRIPT_DIR}/submit_msa_array.sh"
if [ ! -f "$ARRAY_SCRIPT" ]; then
    echo -e "${RED}Error: Array submission script not found: $ARRAY_SCRIPT${NC}"
    exit 1
fi

# Create logs directory if it doesn't exist
mkdir -p logs

# Count total jobs (excluding header)
TOTAL_JOBS=$(tail -n +2 "$CSV_FILE" | wc -l)
echo "Total jobs to submit: $TOTAL_JOBS"

if [ $TOTAL_JOBS -eq 0 ]; then
    echo -e "${RED}Error: No jobs found in CSV file (excluding header)${NC}"
    exit 1
fi

# Calculate split point (half of jobs)
SPLIT_POINT=$(( (TOTAL_JOBS + 1) / 2 ))
PART1_JOBS=$SPLIT_POINT
PART2_JOBS=$(( TOTAL_JOBS - SPLIT_POINT ))

echo "Will split jobs:"
echo "  Part 1 (hns partition): $PART1_JOBS jobs"
echo "  Part 2 (normal partition): $PART2_JOBS jobs"
echo

# Create temporary files for each partition
PART1_FILE="${AF3_BASE_DIR}/msa_array_jobs_part1.tmp"
PART2_FILE="${AF3_BASE_DIR}/msa_array_jobs_part2.tmp"

# Extract header
HEADER=$(head -n 1 "$CSV_FILE")

# Create part 1 file (first half)
echo "$HEADER" > "$PART1_FILE"
tail -n +2 "$CSV_FILE" | head -n $PART1_JOBS >> "$PART1_FILE"

# Create part 2 file (second half)
echo "$HEADER" > "$PART2_FILE"
tail -n +2 "$CSV_FILE" | tail -n $PART2_JOBS >> "$PART2_FILE"

echo -e "${GREEN}Created partition files:${NC}"
echo "  Part 1: $PART1_FILE"
echo "  Part 2: $PART2_FILE"
echo

# Function to validate jobs
validate_jobs() {
    local csv_file=$1
    local part_name=$2
    local invalid_count=0
    local line_num=2  # Start after header
    
    echo -e "${YELLOW}Validating jobs for $part_name...${NC}"
    
    tail -n +2 "$csv_file" | while IFS= read -r job_name; do
        job_name=$(echo "$job_name" | tr -d '\r')
        if [ -z "$job_name" ]; then
            continue
        fi
        
        # Check if job directory exists
        found=false
        if [ -d "${AF3_BASE_DIR}/jobs/${job_name}" ]; then
            if [ -f "${AF3_BASE_DIR}/jobs/${job_name}/alphafold_input.json" ]; then
                found=true
            else
                echo -e "${RED}  Warning: Job $job_name (line $line_num) missing alphafold_input.json${NC}"
                ((invalid_count++))
            fi
        elif [ -d "${AF3_BASE_DIR}/jobs/human_test_set/${job_name}" ]; then
            if [ -f "${AF3_BASE_DIR}/jobs/human_test_set/${job_name}/alphafold_input.json" ]; then
                found=true
            else
                echo -e "${RED}  Warning: Job $job_name (line $line_num) missing alphafold_input.json${NC}"
                ((invalid_count++))
            fi
        else
            echo -e "${RED}  Warning: Job directory not found for $job_name (line $line_num)${NC}"
            ((invalid_count++))
        fi
        
        ((line_num++))
    done
    
    if [ $invalid_count -gt 0 ]; then
        echo -e "${YELLOW}  Found $invalid_count jobs with issues${NC}"
        read -p "Continue with submission? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            return 1
        fi
    else
        echo -e "${GREEN}  All jobs validated successfully${NC}"
    fi
    
    return 0
}

# Function to submit arrays for a partition
submit_partition() {
    local csv_file=$1
    local partition=$2
    local part_name=$3
    local part_jobs=$4
    
    echo -e "${BLUE}Processing $part_name (partition: $partition)${NC}"
    echo "========================================="
    
    # Validate jobs
    if ! validate_jobs "$csv_file" "$part_name"; then
        echo -e "${RED}Skipping $part_name due to validation failures${NC}"
        return 1
    fi
    
    # Calculate number of arrays needed
    local num_arrays=$(( (part_jobs + MAX_ARRAY_SIZE - 1) / MAX_ARRAY_SIZE ))
    echo "Will submit $num_arrays array job(s) for $part_name"
    
    local submitted=0
    for i in $(seq 1 $num_arrays); do
        # Calculate start and end indices for this array
        local start_idx=$(( (i - 1) * MAX_ARRAY_SIZE + 1 ))
        local end_idx=$(( i * MAX_ARRAY_SIZE ))
        if [ $end_idx -gt $part_jobs ]; then
            end_idx=$part_jobs
        fi
        
        # Calculate array size
        local array_size=$(( end_idx - start_idx + 1 ))
        
        echo -e "\n${GREEN}Submitting array $i for $part_name:${NC}"
        echo "  Jobs: $start_idx-$end_idx ($array_size jobs)"
        echo "  Partition: $partition"
        
        # Export environment variables
        export JOB_LIST_FILE="$csv_file"
        export ARRAY_OFFSET=0  # No offset needed since we're using separate files
        
        # Submit with partition specified
        JOB_ID=$(sbatch --parsable --partition=$partition --array=1-${array_size} "$ARRAY_SCRIPT")
        
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}  Successfully submitted with Job ID: $JOB_ID${NC}"
            submitted=$(( submitted + array_size ))
        else
            echo -e "${RED}  Failed to submit array${NC}"
        fi
    done
    
    echo -e "${GREEN}Submitted $submitted jobs for $part_name${NC}\n"
    return 0
}

# Submit both partitions
TOTAL_SUBMITTED=0

# Submit Part 1 to hns partition
if submit_partition "$PART1_FILE" "hns" "Part 1" "$PART1_JOBS"; then
    TOTAL_SUBMITTED=$(( TOTAL_SUBMITTED + PART1_JOBS ))
fi

# Submit Part 2 to normal partition
if submit_partition "$PART2_FILE" "normal" "Part 2" "$PART2_JOBS"; then
    TOTAL_SUBMITTED=$(( TOTAL_SUBMITTED + PART2_JOBS ))
fi

echo "====================================================="
echo -e "${GREEN}Submission complete!${NC}"
echo "Submitted $TOTAL_SUBMITTED out of $TOTAL_JOBS jobs"
echo
echo "Partition files created:"
echo "  $PART1_FILE (hns partition)"
echo "  $PART2_FILE (normal partition)"
echo
echo -e "${YELLOW}Note: These .tmp files must remain accessible during job execution!${NC}"
echo
echo "Monitor jobs with: squeue -u $USER"
echo "Check logs in: logs/"
echo
echo "To clean up .tmp files after all jobs complete:"
echo "  rm $PART1_FILE $PART2_FILE"
