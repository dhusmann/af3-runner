#!/bin/bash

# rclone_to_gdrive.sh - Submit rclone sync job to SLURM
# This script submits an sbatch job to sync organized outputs to Google Drive

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Configuration
BASE_DIR="/scratch/groups/ogozani/alphafold3"
OUTPUT_DIR="${BASE_DIR}/output"
GDRIVE_REMOTE="gozani_labshare_alphafold"
GDRIVE_PATH="alphafold3/output/"
LOG_DIR="${BASE_DIR}/logs"
RCLONE_RATE_FLAGS="--drive-pacer-min-sleep 1s --drive-pacer-burst 4 --tpslimit 4 --tpslimit-burst 4 --transfers 2 --checkers 2"
SBATCH_EXTRA_ARGS="${SBATCH_EXTRA_ARGS:-}"

# Check if output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo -e "${RED}Error: Output directory not found: $OUTPUT_DIR${NC}"
    echo "Run sync_organize_outputs.sh first"
    exit 1
fi

# Create logs directory if needed
mkdir -p "$LOG_DIR"

# Display sync information (without slow file counting)
echo "Preparing rclone sync job..."

echo "=== OUTPUT SYNC ==="
echo "Source: $OUTPUT_DIR → ${GDRIVE_REMOTE}:${GDRIVE_PATH}"
echo

echo "=== SCRIPTS & ANALYSIS BACKUP ==="  
echo "Scripts: $BASE_DIR → ${GDRIVE_REMOTE}:alphafold3/scripts/"
echo "Tools: $BASE_DIR/tools → ${GDRIVE_REMOTE}:alphafold3/scripts/tools/"
echo "Analysis: $BASE_DIR/analysis → ${GDRIVE_REMOTE}:alphafold3/analysis/"
echo

# Create the sbatch script
SBATCH_SCRIPT="${BASE_DIR}/rclone_sync_job.sh"
cat > "$SBATCH_SCRIPT" << 'EOF'
#!/bin/bash
#SBATCH --job-name=SyncAF3_Gdrive
#SBATCH --output=logs/Rclone_%j.out
#SBATCH --error=logs/Rclone_%j.err
#SBATCH --partition=hns
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=dhusmann@stanford.edu

# Load modules
ml system
ml rclone/1.59.1

# Set variables from environment
OUTPUT_DIR="${SYNC_OUTPUT_DIR}"
GDRIVE_REMOTE="${SYNC_GDRIVE_REMOTE}"
GDRIVE_PATH="${SYNC_GDRIVE_PATH}"
BASE_DIR="${SYNC_BASE_DIR}"

echo "Starting rclone sync at $(date)"
echo "Base directory: $BASE_DIR"
echo

# Step 1: Backup scripts (top-level files, excluding logs)
echo "=== STEP 1: Backing up scripts and configuration files ==="
echo "Source: $BASE_DIR (top-level files)"
echo "Destination: ${GDRIVE_REMOTE}:alphafold3/scripts/"
rclone copy -v "$BASE_DIR/" "${GDRIVE_REMOTE}:alphafold3/scripts/" \
    --exclude "logs/**" \
    --exclude "*.out" \
    --exclude "*.err" \
    --exclude "*.log" \
    --exclude "*.tmp" \
    --exclude "jobs/**" \
    --exclude "output/**" \
    --exclude "analysis/**" \
    --exclude "tools/**" \
    --exclude "rclone_sync_job.sh" \
    --max-depth 1 \
    $RCLONE_RATE_FLAGS \
    2>&1 | tee -a rclone_output.tmp

echo
echo "=== STEP 2: Backing up tools directory ==="
echo "Source: $BASE_DIR/tools/"
echo "Destination: ${GDRIVE_REMOTE}:alphafold3/scripts/tools/"
if [ -d "$BASE_DIR/tools" ]; then
    rclone copy -v "$BASE_DIR/tools/" "${GDRIVE_REMOTE}:alphafold3/scripts/tools/" \
        $RCLONE_RATE_FLAGS \
        2>&1 | tee -a rclone_output.tmp
else
    echo "Tools directory not found, skipping..."
fi

echo
echo "=== STEP 3: Backing up analysis directory ==="
echo "Source: $BASE_DIR/analysis/"
echo "Destination: ${GDRIVE_REMOTE}:alphafold3/analysis/"
if [ -d "$BASE_DIR/analysis" ]; then
    rclone copy -v "$BASE_DIR/analysis/" "${GDRIVE_REMOTE}:alphafold3/analysis/" \
        $RCLONE_RATE_FLAGS \
        2>&1 | tee -a rclone_output.tmp
else
    echo "Analysis directory not found, skipping..."
fi

echo
echo "=== STEP 4: Syncing output data ==="
echo "Source: $OUTPUT_DIR"
echo "Destination: ${GDRIVE_REMOTE}:${GDRIVE_PATH}"

# Run rclone with quota detection
# Use --filter flags for deterministic rule processing
# Include MSA archives and seed archives, exclude uncompressed seed directories
rclone copy -v "$OUTPUT_DIR/" "${GDRIVE_REMOTE}:${GDRIVE_PATH}" \
    --filter "+ msa/**" \
    --filter "+ **/seeds.tar.gz" \
    --filter "+ **/seeds.tar" \
    --filter "- **/seed-*/**" \
    --filter "- **/seed-*" \
    --filter "+ **" \
    $RCLONE_RATE_FLAGS \
    2>&1 | tee -a rclone_output.tmp

# Check for quota errors
if grep -q "teamDriveFileLimitExceeded\|userRateLimitExceeded" rclone_output.tmp; then
    echo
    echo "ERROR: Google Drive quota exceeded!"
    echo "Sync stopped at $(date)"
    
    # Send email notification
    mail -s "AlphaFold3 Sync: Google Drive Quota Exceeded" dhusmann@stanford.edu << MAIL_END
The AlphaFold3 output sync to Google Drive has been stopped due to quota limits.

Job ID: $SLURM_JOB_ID
Time: $(date)
Files synced before quota: $(grep -c "Copied (new)" rclone_output.tmp || echo "0")

The sync job will be automatically rescheduled to run in 24 hours.

To check sync status:
  cd $BASE_DIR
  ./check_sync_status.sh

To manually restart sync later:
  ./rclone_to_gdrive.sh
MAIL_END
    
    # Schedule a new job for 24 hours later
    echo "Scheduling retry in 24 hours..."
    
    # Create retry script
    cat > "${BASE_DIR}/rclone_retry.sh" << RETRY_SCRIPT
#!/bin/bash
cd "$BASE_DIR"
./rclone_to_gdrive.sh
RETRY_SCRIPT
    
    chmod +x "${BASE_DIR}/rclone_retry.sh"
    
    # Submit with 24 hour delay
    sbatch $SBATCH_EXTRA_ARGS --begin=now+24hours --job-name=SyncAF3_Retry "${BASE_DIR}/rclone_retry.sh"
    
    rm rclone_output.tmp
    exit 1
else
    # Success
    echo
    echo "Sync completed successfully at $(date)"
    
    # Count what was synced for each step
    SCRIPTS_COPIED=$(grep -A 20 "STEP 1:" rclone_output.tmp | grep -c "Copied (new)" || echo "0")
    TOOLS_COPIED=$(grep -A 20 "STEP 2:" rclone_output.tmp | grep -c "Copied (new)" || echo "0") 
    ANALYSIS_COPIED=$(grep -A 20 "STEP 3:" rclone_output.tmp | grep -c "Copied (new)" || echo "0")
    OUTPUT_COPIED=$(grep -A 1000 "STEP 4:" rclone_output.tmp | grep -c "Copied (new)" || echo "0")
    TOTAL_COPIED=$(grep -c "Copied (new)" rclone_output.tmp || echo "0")
    
    echo "=== SYNC SUMMARY ==="
    echo "Scripts backed up: $SCRIPTS_COPIED files"
    echo "Tools backed up: $TOOLS_COPIED files" 
    echo "Analysis backed up: $ANALYSIS_COPIED files"
    echo "Output files synced: $OUTPUT_COPIED files"
    echo "Total files processed: $TOTAL_COPIED files"
    
    rm rclone_output.tmp
fi
EOF

# Make script executable
chmod +x "$SBATCH_SCRIPT"

# Export variables for the sbatch job
export SYNC_OUTPUT_DIR="$OUTPUT_DIR"
export SYNC_GDRIVE_REMOTE="$GDRIVE_REMOTE"
export SYNC_GDRIVE_PATH="$GDRIVE_PATH"
export SYNC_BASE_DIR="$BASE_DIR"
export SBATCH_EXTRA_ARGS

# Submit the job
echo -e "${GREEN}Submitting rclone sync job to SLURM...${NC}"
JOB_ID=$(sbatch $SBATCH_EXTRA_ARGS --export=ALL "$SBATCH_SCRIPT" | awk '{print $4}')

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Successfully submitted sync job with ID: $JOB_ID${NC}"
    echo
    echo "Monitor progress with:"
    echo "  squeue -j $JOB_ID"
    echo "  tail -f logs/Rclone_${JOB_ID}.out"
    echo
    echo "The job will:"
    echo "  - Backup scripts to Google Drive (alphafold3/scripts/)"
    echo "  - Backup tools to Google Drive (alphafold3/scripts/tools/)"
    echo "  - Backup analysis to Google Drive (alphafold3/analysis/)"
    echo "  - Sync all output files to Google Drive (alphafold3/output/)"
    echo "  - Stop if quota is exceeded and email you"
    echo "  - Automatically retry in 24 hours if quota hit"
else
    echo -e "${RED}Failed to submit job${NC}"
    exit 1
fi
