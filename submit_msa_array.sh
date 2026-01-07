#!/bin/bash
#SBATCH --job-name="af3_msa_array"
#SBATCH --output=logs/%A_%a_MSA_array_v2.out
#SBATCH --error=logs/%A_%a_MSA_array_v2.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=72GB
#SBATCH --time=8:00:00

# Record start time
START_TIME=$(date +%s)
START_TIME_HUMAN=$(date "+%Y-%m-%d %H:%M:%S")

# Base directory for AlphaFold 3 resources
export AF3_BASE_DIR=/scratch/groups/ogozani/alphafold3
LOG_FILE="${AF3_BASE_DIR}/logged_folding_jobs.csv"

# Get the job list file from environment variable (set by wrapper script)
JOB_LIST_FILE="${JOB_LIST_FILE:-${AF3_BASE_DIR}/msa_array_jobs.csv}"

# Calculate the actual line number (accounting for header)
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))

# Get the job name from the CSV file
JOB_NAME=$(sed -n "${LINE_NUM}p" "$JOB_LIST_FILE" | tr -d '\r')

if [ -z "$JOB_NAME" ]; then
    echo "Error: Could not find job at line $LINE_NUM in $JOB_LIST_FILE"
    exit 1
fi

# Print job info
echo "Running AlphaFold 3 MSA generation for folder: ${JOB_NAME}"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Start time: ${START_TIME_HUMAN}"

# Check if job exists in jobs/ or jobs/human_test_set/
if [ -d "${AF3_BASE_DIR}/jobs/${JOB_NAME}" ]; then
    JOB_DIR="${AF3_BASE_DIR}/jobs/${JOB_NAME}"
elif [ -d "${AF3_BASE_DIR}/jobs/human_test_set/${JOB_NAME}" ]; then
    JOB_DIR="${AF3_BASE_DIR}/jobs/human_test_set/${JOB_NAME}"
else
    echo "Error: Job directory not found for ${JOB_NAME}"
    echo "Checked: ${AF3_BASE_DIR}/jobs/${JOB_NAME}"
    echo "Checked: ${AF3_BASE_DIR}/jobs/human_test_set/${JOB_NAME}"
    exit 1
fi

# Check if alphafold_input.json exists
if [ ! -f "${JOB_DIR}/alphafold_input.json" ]; then
    echo "Error: alphafold_input.json not found in ${JOB_DIR}"
    exit 1
fi

echo "Found job directory: ${JOB_DIR}"

# Set up all the paths
export AF3_IMAGE=${AF3_BASE_DIR}/alphafold3_resources/image/alphafold3.sif
export AF3_INPUT_DIR=${JOB_DIR}
export AF3_OUTPUT_DIR=${JOB_DIR}/output_msa
export AF3_DATABASES_DIR=${AF3_BASE_DIR}/alphafold3_resources/databases

# Create output directory if it doesn't exist
mkdir -p ${AF3_OUTPUT_DIR}

# Function to update the log file
update_log_file() {
    local exit_code=$1
    local end_time=$(date +%s)
    local end_time_human=$(date "+%Y-%m-%d %H:%M:%S")
    local duration=$((end_time - START_TIME))
    
    # Convert duration to human-readable format
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))
    local duration_human=$(printf "%02d:%02d:%02d" $hours $minutes $seconds)
    
    # Update log file (implementation same as before)
    local lock_file="${AF3_BASE_DIR}/.log_update.lock"
    local max_wait=30
    local wait_count=0
    
    while [ -f "$lock_file" ] && [ $wait_count -lt $max_wait ]; do
        sleep 1
        ((wait_count++))
    done
    
    touch "$lock_file"
    
    # Add MSA completion info to log
    if [ -f "$LOG_FILE" ]; then
        local temp_log="${AF3_BASE_DIR}/.logged_folding_jobs_${SLURM_JOB_ID}.tmp"
        awk -F',' -v job_name="${JOB_NAME}" -v job_id="${SLURM_JOB_ID}" -v duration="${duration_human}" \
            -v exit_code="${exit_code}" -v end_time="${end_time_human}" \
            'BEGIN {OFS=","} 
            $1 == job_name && $2 == job_id && $3 == "MSA" {
                # Update MSA stage completion
                $5 = end_time
                $6 = duration
                $7 = exit_code
            }
            {print}' "$LOG_FILE" > "$temp_log"
        mv "$temp_log" "$LOG_FILE"
    fi
    
    rm -f "$lock_file"
    
    echo "MSA stage completion logged: Duration=${duration_human}, Exit code=${exit_code}"
}

# Set trap to update log file on exit
trap 'update_log_file $?' EXIT

# Run AlphaFold 3 - MSA generation only
singularity exec \
     --bind $AF3_INPUT_DIR:/root/af_input \
     --bind $AF3_OUTPUT_DIR:/root/af_output \
     --bind $AF3_DATABASES_DIR:/root/public_databases \
     $AF3_IMAGE \
     python /app/alphafold3/run_alphafold.py \
     --json_path=/root/af_input/alphafold_input.json \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output \
     --norun_inference
