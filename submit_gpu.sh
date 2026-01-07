#!/bin/bash
#SBATCH --job-name="af3_gpu"
#SBATCH --output=logs/%A_GPU-only.out
#SBATCH --error=logs/%A_GPU-only.err
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --gpus=1
#SBATCH -C "GPU_SKU:L40S|GPU_SKU:H100_SXM5"
#SBATCH --time=2:00:00

# Note: The partition is now set via command line in submit_dist.sh using -p flag

# Record start time
START_TIME=$(date +%s)
START_TIME_HUMAN=$(date "+%Y-%m-%d %H:%M:%S")

# Get the folder name from command line argument
FOLDER_NAME=${1}

# Print job info
echo "Running AlphaFold 3 GPU inference for folder: ${FOLDER_NAME}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo "Start time: ${START_TIME_HUMAN}"

# Base directory for AlphaFold 3 resources
export AF3_BASE_DIR=/scratch/groups/ogozani/alphafold3
LOG_FILE="${AF3_BASE_DIR}/logged_folding_jobs.csv"

# Function to find job directory
find_job_directory() {
    local job_name=$1
    
    # Check primary location
    if [ -d "${AF3_BASE_DIR}/jobs/$job_name" ]; then
        echo "${AF3_BASE_DIR}/jobs/$job_name"
        return 0
    fi
    
    # Check human_test_set location
    if [ -d "${AF3_BASE_DIR}/jobs/human_test_set/$job_name" ]; then
        echo "${AF3_BASE_DIR}/jobs/human_test_set/$job_name"
        return 0
    fi
    
    # Not found
    echo "Error: Job directory not found for $job_name" >&2
    return 1
}

# Find the job directory
JOB_DIR=$(find_job_directory "${FOLDER_NAME}")
if [ $? -ne 0 ]; then
    echo "Fatal error: Cannot find job directory for ${FOLDER_NAME}"
    exit 1
fi

echo "Using job directory: $JOB_DIR"

# Set up all the paths
export AF3_IMAGE=${AF3_BASE_DIR}/alphafold3_resources/image/alphafold3.sif
export AF3_MSA_OUTPUT_DIR=${JOB_DIR}/output_msa
export AF3_OUTPUT_DIR=${JOB_DIR}/output
export AF3_MODEL_PARAMETERS_DIR=${AF3_BASE_DIR}/alphafold3_resources/weights
export AF3_DATABASES_DIR=${AF3_BASE_DIR}/alphafold3_resources/databases

# Create final output directory
mkdir -p ${AF3_OUTPUT_DIR}

# Detect if this is a single-protein + ligand job by inspecting alphafold_input.json
IS_SINGLE=$(singularity exec \
     --nv \
     --bind ${JOB_DIR}:/job \
     $AF3_IMAGE \
     python - <<'PY'
import json, sys
try:
    d=json.load(open('/job/alphafold_input.json'))
    nprot=sum(1 for e in d.get('sequences',[]) if 'protein' in e)
    nlig=sum(1 for e in d.get('sequences',[]) if 'ligand' in e)
    print(1 if (nprot==1 and nlig==1) else 0)
except Exception:
    print(0)
PY)

if [ "$IS_SINGLE" = "1" ]; then
    echo "Single-protein+ligand detected → running full data pipeline."
    JSON_PATH="/root/af_job/alphafold_input.json"
else
    # Find the augmented JSON file from MSA stage (prefer consolidated file)
    if [ -f "${AF3_MSA_OUTPUT_DIR}/alphafold_input_with_msa.json" ]; then
        JSON_PATH="/root/af_msa_output/alphafold_input_with_msa.json"
    else
        JSON_BASENAME=$(find ${AF3_MSA_OUTPUT_DIR} -maxdepth 1 -name "*.json" -type f | head -1 | xargs -I{} basename {})
        JSON_PATH="/root/af_msa_output/${JSON_BASENAME}"
    fi
    if [ -z "$JSON_PATH" ]; then
        echo "Error: No augmented JSON file found in ${AF3_MSA_OUTPUT_DIR}"
        exit 1
    fi
    echo "Using augmented JSON: ${JSON_PATH}"
fi

# Function to update the log file
update_log_file() {
    local exit_code=$1
    local end_time=$(date +%s)
    local end_time_human=$(date "+%Y-%m-%d %H:%M:%S")
    local duration=$((end_time - START_TIME))
    
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))
    local duration_human=$(printf "%02d:%02d:%02d" $hours $minutes $seconds)
    
    local lock_file="${AF3_BASE_DIR}/.log_update.lock"
    local max_wait=30
    local wait_count=0
    
    while [ -f "$lock_file" ] && [ $wait_count -lt $max_wait ]; do
        sleep 1
        ((wait_count++))
    done
    
    touch "$lock_file"
    
    if [ -f "$LOG_FILE" ]; then
        local temp_log="${AF3_BASE_DIR}/.logged_folding_jobs_${SLURM_JOB_ID}.tmp"
        awk -F',' -v job_name="${FOLDER_NAME}" -v job_id="${SLURM_JOB_ID}" -v duration="${duration_human}" \
            -v exit_code="${exit_code}" -v end_time="${end_time_human}" \
            'BEGIN {OFS=","} 
            $1 == job_name && $2 == job_id && $3 == "GPU" {
                # Update GPU stage completion
                $5 = end_time
                $6 = duration
                $7 = exit_code
            }
            {print}' "$LOG_FILE" > "$temp_log"
        mv "$temp_log" "$LOG_FILE"
    fi
    
    rm -f "$lock_file"
    
    echo "GPU stage completion logged: Duration=${duration_human}, Exit code=${exit_code}"
}

trap 'update_log_file $?' EXIT

# Tune CPU threading for data pipeline when cpus-per-task is provided
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Run AlphaFold 3
singularity exec \
     --nv \
     $( [ "$IS_SINGLE" = "1" ] && echo "--bind ${JOB_DIR}:/root/af_job" ) \
     --bind ${AF3_MSA_OUTPUT_DIR}:/root/af_msa_output \
     --bind $AF3_OUTPUT_DIR:/root/af_output \
     --bind $AF3_MODEL_PARAMETERS_DIR:/root/models \
     --bind $AF3_DATABASES_DIR:/root/public_databases \
     $AF3_IMAGE \
     python /app/alphafold3/run_alphafold.py \
     --json_path=${JSON_PATH} \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output \
     $( [ "$IS_SINGLE" = "1" ] && echo "" || echo "--norun_data_pipeline" )
