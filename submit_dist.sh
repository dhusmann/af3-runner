#!/bin/bash

# Default values
MAX_GPU_JOBS=50
MAX_HNS_JOBS=50
CSV_FILE="folding_jobs.csv"
LOG_FILE="logged_folding_jobs.csv"
TEMP_CSV="folding_jobs.tmp"
WAITING_CSV="waiting_for_msa.csv"
DEBUG_MODE=0

# Function to display usage
usage() {
    echo "Usage: $0 [-g MAX_GPU] [-n MAX_HNS] [-d]"
    echo "  -g MAX_GPU     Maximum number of GPU partition jobs (default: 50)"
    echo "  -n MAX_HNS     Maximum number of HNS partition jobs (default: 50)"
    echo "  -d             Enable debug mode"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "g:n:dh" opt; do
    case $opt in
        g)
            MAX_GPU_JOBS=$OPTARG
            ;;
        n)
            MAX_HNS_JOBS=$OPTARG
            ;;
        d)
            DEBUG_MODE=1
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

# Function to check if MSA arrays are still running
check_msa_arrays_running() {
    local running_count=$(squeue -u $USER -n af3_msa_array -h | wc -l)
    return $running_count
}

# Run batch_reuse_msa.py to update MSA status and identify waiting jobs
echo "=== Running MSA reuse analysis ==="

# Check if MSA array jobs are still running
check_msa_arrays_running
MSA_RUNNING=$?

if [ $MSA_RUNNING -gt 0 ]; then
    echo "Note: $MSA_RUNNING MSA array jobs are still running"
    echo "Some jobs in waiting_for_msa.csv may become ready as these complete"
fi

# Load python module
ml python/3.9.0 2>/dev/null || module load python/3.9.0 2>/dev/null || true

# Run batch_reuse_msa.py
if [ -f "batch_reuse_msa.py" ]; then
    python batch_reuse_msa.py
    if [ $? -ne 0 ]; then
        echo "Warning: batch_reuse_msa.py failed, continuing with existing files"
    fi
else
    echo "Warning: batch_reuse_msa.py not found, skipping MSA reuse"
fi

echo ""

# Check if CSV file exists
if [ ! -f "$CSV_FILE" ]; then
    echo "Error: $CSV_FILE not found!"
    exit 1
fi

# Create log file header if it doesn't exist
if [ ! -f "$LOG_FILE" ]; then
    echo "job_name,slurm_job_id,stage,submission_date,completion_date,duration,exit_code" > "$LOG_FILE"
fi

# Function to count jobs by partition
count_partition_jobs() {
    local partition=$1
    squeue -u $USER -h -p $partition | grep "af3_gpu" | wc -l
}

# Get current job counts for both partitions
CURRENT_GPU_JOBS=$(count_partition_jobs "gpu")
CURRENT_HNS_JOBS=$(count_partition_jobs "hns")

echo "Current GPU partition jobs: $CURRENT_GPU_JOBS / $MAX_GPU_JOBS"
echo "Current HNS partition jobs: $CURRENT_HNS_JOBS / $MAX_HNS_JOBS"

# Calculate available slots for each partition
AVAILABLE_GPU_SLOTS=$((MAX_GPU_JOBS - CURRENT_GPU_JOBS))
AVAILABLE_HNS_SLOTS=$((MAX_HNS_JOBS - CURRENT_HNS_JOBS))

if [ $AVAILABLE_GPU_SLOTS -lt 0 ]; then
    AVAILABLE_GPU_SLOTS=0
fi
if [ $AVAILABLE_HNS_SLOTS -lt 0 ]; then
    AVAILABLE_HNS_SLOTS=0
fi

# Store original available slots for tracking
ORIGINAL_GPU_SLOTS=$AVAILABLE_GPU_SLOTS
ORIGINAL_HNS_SLOTS=$AVAILABLE_HNS_SLOTS
TOTAL_AVAILABLE=$((AVAILABLE_GPU_SLOTS + AVAILABLE_HNS_SLOTS))

if [ $TOTAL_AVAILABLE -le 0 ]; then
    echo "Already at or above job limits for both partitions. No jobs will be submitted."
    exit 0
fi

echo "Available GPU partition slots: $AVAILABLE_GPU_SLOTS"
echo "Available HNS partition slots: $AVAILABLE_HNS_SLOTS"
echo "Total available slots: $TOTAL_AVAILABLE"

# Function to find job directory
find_job_directory() {
    local job_name=$1
    
    if [ -d "jobs/$job_name" ]; then
        echo "jobs/$job_name"
        return 0
    fi
    
    if [ -d "jobs/human_test_set/$job_name" ]; then
        echo "jobs/human_test_set/$job_name"
        return 0
    fi
    
    return 1
}

# Function to check if MSA is complete
msa_is_complete() {
    local job_name=$1
    local job_dir=$(find_job_directory "$job_name")
    
    if [ -z "$job_dir" ]; then
        return 1
    fi
    
    local msa_output_dir="${job_dir}/output_msa"
    
    if [ -d "$msa_output_dir" ] && [ -n "$(find "$msa_output_dir" -name "*.json" -type f 2>/dev/null | head -1)" ]; then
        return 0
    else
        return 1
    fi
}

# Function to check if GPU is complete
gpu_is_complete() {
    local job_name=$1
    local job_dir=$(find_job_directory "$job_name")
    
    if [ -z "$job_dir" ]; then
        return 1
    fi
    
    local gpu_output_dir="${job_dir}/output"
    
    if [ -d "$gpu_output_dir" ] && [ -n "$(ls -A "$gpu_output_dir" 2>/dev/null)" ]; then
        return 0
    else
        return 1
    fi
}

# Function to detect single-protein-ligand job (1 protein, 1 ligand)
is_single_protein_ligand() {
    local job_name=$1
    local job_dir=$(find_job_directory "$job_name")
    if [ -z "$job_dir" ]; then
        echo 0
        return 0
    fi
    python - << 'PY'
import json, sys, os
job_dir = os.environ.get('JOB_DIR')
try:
    with open(os.path.join(job_dir, 'alphafold_input.json')) as fh:
        d = json.load(fh)
except Exception:
    print(0)
    sys.exit(0)
nprot = sum(1 for e in d.get('sequences', []) if 'protein' in e)
nlig  = sum(1 for e in d.get('sequences', []) if 'ligand' in e)
print(1 if (nprot==1 and nlig==1) else 0)
PY
}

# Function to determine which partition to use
get_next_partition() {
    # Check actual remaining slots (not the decremented ones)
    local gpu_remaining=$((ORIGINAL_GPU_SLOTS - SUBMITTED_GPU))
    local hns_remaining=$((ORIGINAL_HNS_SLOTS - SUBMITTED_HNS))
    
    # If we've hit the total limit, stop
    if [ $SUBMITTED -ge $TOTAL_AVAILABLE ]; then
        echo "none"
        return
    fi
    
    # Prefer the partition with more remaining slots
    if [ $gpu_remaining -gt $hns_remaining ]; then
        if [ $gpu_remaining -gt 0 ]; then
            echo "gpu"
        elif [ $hns_remaining -gt 0 ]; then
            echo "hns"
        else
            echo "none"
        fi
    else
        if [ $hns_remaining -gt 0 ]; then
            echo "hns"
        elif [ $gpu_remaining -gt 0 ]; then
            echo "gpu"
        else
            echo "none"
        fi
    fi
}

# Function to process a single job
process_job() {
    local folder_name=$1
    local from_waiting=$2
    
    if [ $DEBUG_MODE -eq 1 ]; then
        echo "Processing: $folder_name (from_waiting=$from_waiting)"
    fi
    
    # Check if we've already submitted enough jobs
    if [ $SUBMITTED -ge $TOTAL_AVAILABLE ]; then
        return 3  # Stop processing
    fi
    
    # Find job directory
    job_dir=$(find_job_directory "$folder_name")
    
    if [ -z "$job_dir" ]; then
        return 1  # Keep in CSV
    fi
    
    # Check if alphafold_input.json exists
    if [ ! -f "${job_dir}/alphafold_input.json" ]; then
        return 1  # Keep in CSV
    fi
    
    # Check if GPU is complete
    if gpu_is_complete "$folder_name"; then
        ((REMOVED_COMPLETE++))
        if [ $DEBUG_MODE -eq 1 ]; then
            echo "  Removing (complete): $folder_name"
        fi
        return 2  # Remove from CSV
    fi
    
    # Detect single-protein-ligand jobs which do not require precomputed MSA
    export JOB_DIR
    JOB_DIR=$(find_job_directory "$folder_name")
    SINGLE=$(JOB_DIR="$JOB_DIR" is_single_protein_ligand "$folder_name")
    
    # Check if MSA is complete for multi-chain jobs only
    if [ "$SINGLE" != "1" ]; then
        if ! msa_is_complete "$folder_name"; then
            ((SKIPPED_NO_MSA++))
            return 1  # Keep in CSV
        fi
    fi
    
    # Determine which partition to use
    PARTITION=$(get_next_partition)
    
    if [ "$PARTITION" = "none" ]; then
        # No slots available
        return 3  # Stop processing
    fi
    
    # Submit the job to the selected partition
    echo "Submitting GPU job for: $folder_name (partition: $PARTITION)"
    # For single-protein-ligand jobs, allocate more CPUs and time since we run data pipeline
    if [ "$SINGLE" = "1" ]; then
        SBATCH_OUTPUT=$(sbatch -p "$PARTITION" -c 12 --time=08:00:00 submit_gpu.sh "$folder_name" 2>&1)
    else
        SBATCH_OUTPUT=$(sbatch -p "$PARTITION" submit_gpu.sh "$folder_name" 2>&1)
    fi
    
    if [ $? -eq 0 ]; then
        JOB_ID=$(echo "$SBATCH_OUTPUT" | grep -oE '[0-9]+$')
        TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
        echo "${folder_name},${JOB_ID},GPU,${TIMESTAMP},,," >> "$LOG_FILE"
        ((SUBMITTED++))
        
        if [ "$PARTITION" = "gpu" ]; then
            ((SUBMITTED_GPU++))
        else
            ((SUBMITTED_HNS++))
        fi
        
        if [ "$from_waiting" -eq 1 ]; then
            ((SUBMITTED_FROM_WAITING++))
        fi
        
        # Check if we've hit our limit
        if [ $SUBMITTED -ge $TOTAL_AVAILABLE ]; then
            echo "Reached submission limit ($SUBMITTED jobs submitted)"
            return 3  # Stop processing
        fi
        
        # Remove from CSV since successfully submitted
        return 2
    else
        echo "  Error: $SBATCH_OUTPUT"
        return 1  # Keep in CSV to retry later
    fi
}

# Initialize counters
SUBMITTED=0
SUBMITTED_GPU=0
SUBMITTED_HNS=0
SUBMITTED_FROM_WAITING=0
SKIPPED_NO_MSA=0
REMOVED_COMPLETE=0
LINES_PROCESSED=0
STOP_PROCESSING=0

# Create temporary files
> "$TEMP_CSV"
TEMP_WAITING="waiting_for_msa.tmp"
> "$TEMP_WAITING"

echo ""
echo "Processing jobs (will stop after submitting $TOTAL_AVAILABLE jobs)..."

# First, process waiting_for_msa.csv if it exists
if [ -f "$WAITING_CSV" ] && [ $SUBMITTED -lt $TOTAL_AVAILABLE ]; then
    echo ""
    echo "Checking jobs waiting for MSA..."
    
    while IFS= read -r line; do
        # Skip empty lines and header
        if [[ -z "$line" ]] || [[ "$line" =~ ^[[:space:]]*$ ]] || [[ "$line" =~ input_folder_name ]]; then
            if [[ "$line" =~ input_folder_name ]]; then
                echo "$line" >> "$TEMP_WAITING"
            fi
            continue
        fi
        
        # Strip CRs from CRLF lines, remove quotes, trim whitespace
        folder_name=$(printf "%s" "$line" | tr -d '\r"' | xargs)
        
        process_job "$folder_name" 1
        status=$?
        
        if [ $status -eq 3 ]; then
            # Hit limit, keep remaining jobs
            echo "$line" >> "$TEMP_WAITING"
            STOP_PROCESSING=1
            break
        elif [ $status -eq 1 ]; then
            # Keep in waiting list
            echo "$line" >> "$TEMP_WAITING"
        fi
        # status 2 means completed or submitted, don't add to temp file
        
    done < "$WAITING_CSV"
    
    # Replace waiting CSV
    mv "$TEMP_WAITING" "$WAITING_CSV"
fi

# Process main CSV
if [ $STOP_PROCESSING -eq 0 ] && [ $SUBMITTED -lt $TOTAL_AVAILABLE ]; then
    while IFS= read -r line; do
        ((LINES_PROCESSED++))
        
        # Skip empty lines
        if [[ -z "$line" ]] || [[ "$line" =~ ^[[:space:]]*$ ]]; then
            continue
        fi
        
        # Check for header
        if [[ "$line" =~ folder ]] || [[ "$line" =~ input_folder_name ]] || [ $LINES_PROCESSED -eq 1 ]; then
            echo "$line" >> "$TEMP_CSV"
            continue
        fi
        
        # Check if we've submitted enough jobs
        if [ $SUBMITTED -ge $TOTAL_AVAILABLE ]; then
            # Just copy remaining lines without processing
            echo "$line" >> "$TEMP_CSV"
            continue
        fi
        
        # Process this job
        # Strip CRs from CRLF lines, remove quotes, trim whitespace
        folder_name=$(printf "%s" "$line" | tr -d '\r"' | xargs)
        
        process_job "$folder_name" 0
        status=$?
        
        if [ $status -eq 3 ]; then
            # Hit limit
            echo "$line" >> "$TEMP_CSV"
            STOP_PROCESSING=1
        elif [ $status -eq 1 ]; then
            # Keep in CSV
            echo "$line" >> "$TEMP_CSV"
        fi
        # status 2 means completed or submitted, don't add to temp file
        
    done < "$CSV_FILE"
fi

# Replace original CSV
mv "$TEMP_CSV" "$CSV_FILE"

# Summary
echo ""
echo "=== Submission Summary ==="
echo "Lines processed: $LINES_PROCESSED"
echo "Total GPU jobs submitted: $SUBMITTED (Limit: $TOTAL_AVAILABLE)"
if [ $SUBMITTED_GPU -gt 0 ]; then
    echo "  - To GPU partition: $SUBMITTED_GPU / $ORIGINAL_GPU_SLOTS available"
fi
if [ $SUBMITTED_HNS -gt 0 ]; then
    echo "  - To HNS partition: $SUBMITTED_HNS / $ORIGINAL_HNS_SLOTS available"
fi
if [ $SUBMITTED_FROM_WAITING -gt 0 ]; then
    echo "  - From waiting list: $SUBMITTED_FROM_WAITING"
fi
echo "Jobs removed (complete): $REMOVED_COMPLETE"
echo "Jobs skipped (MSA incomplete): $SKIPPED_NO_MSA"

# Recalculate current totals
FINAL_GPU_JOBS=$(count_partition_jobs "gpu")
FINAL_HNS_JOBS=$(count_partition_jobs "hns")
echo ""
echo "Current GPU partition jobs: $FINAL_GPU_JOBS"
echo "Current HNS partition jobs: $FINAL_HNS_JOBS"
echo "Total GPU jobs in queue: $((FINAL_GPU_JOBS + FINAL_HNS_JOBS))"
echo ""

# Quick count of remaining jobs
REMAINING=$(grep -v "folder" "$CSV_FILE" 2>/dev/null | grep -v "^$" | wc -l)
echo "Jobs remaining in main CSV: $REMAINING"

if [ -f "$WAITING_CSV" ]; then
    WAITING_COUNT=$(grep -v "input_folder_name" "$WAITING_CSV" 2>/dev/null | grep -v "^$" | wc -l)
    echo "Jobs waiting for MSA: $WAITING_COUNT"
fi

# Check MSA array status again
check_msa_arrays_running
MSA_STILL_RUNNING=$?
if [ $MSA_STILL_RUNNING -gt 0 ]; then
    echo ""
    echo "Note: $MSA_STILL_RUNNING MSA array jobs still running"
fi
