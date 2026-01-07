#!/bin/bash

# Utility script to monitor MSA array jobs

# Configuration
AF3_BASE_DIR="/scratch/groups/ogozani/alphafold3"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}AlphaFold 3 MSA Array Job Monitor${NC}"
echo "=================================="

# Function to count job states by partition
count_jobs() {
    local state=$1
    local partition=$2
    if [ -z "$partition" ]; then
        squeue -u $USER -n af3_msa_array -t $state -h | wc -l
    else
        squeue -u $USER -n af3_msa_array -t $state -p $partition -h | wc -l
    fi
}

# Show summary
echo -e "\n${GREEN}Job Summary:${NC}"
echo "  Running (total):    $(count_jobs RUNNING)"
echo "    - normal:        $(count_jobs RUNNING normal)"
echo "    - hns:           $(count_jobs RUNNING hns)"
echo "  Pending (total):    $(count_jobs PENDING)"
echo "    - normal:        $(count_jobs PENDING normal)"
echo "    - hns:           $(count_jobs PENDING hns)"
echo "  Completing:         $(count_jobs COMPLETING)"

# Show detailed array job status by partition
echo -e "\n${GREEN}Array Jobs (normal partition):${NC}"
squeue -u $USER -n af3_msa_array -p normal -o "%.18i %.9P %.50j %.8u %.2t %.10M %.6D %R" | grep -E "(JOBID|af3_msa_array\[)" | head -10

echo -e "\n${GREEN}Array Jobs (hns partition):${NC}"
squeue -u $USER -n af3_msa_array -p hns -o "%.18i %.9P %.50j %.8u %.2t %.10M %.6D %R" | grep -E "(JOBID|af3_msa_array\[)" | head -10

# Check for recent failures
echo -e "\n${GREEN}Recent Job Failures (last 24h):${NC}"
FAILURES=$(sacct -u $USER -S $(date -d '24 hours ago' +%Y-%m-%d) -n -X --name=af3_msa_array --state=FAILED,CANCELLED,TIMEOUT --format=JobID,JobName%50,State,ExitCode,Elapsed)

if [ -z "$FAILURES" ]; then
    echo "  No failures found"
else
    echo "$FAILURES"
fi

# Check log files for errors
echo -e "\n${GREEN}Recent Errors in Logs:${NC}"
if [ -d "logs" ]; then
    RECENT_ERRORS=$(find logs -name "*_MSA.err" -mtime -1 -exec grep -l "Error:" {} \; 2>/dev/null | head -10)
    if [ -z "$RECENT_ERRORS" ]; then
        echo "  No error logs found"
    else
        echo "  Error logs found in:"
        echo "$RECENT_ERRORS" | sed 's/^/    /'
    fi
else
    echo "  Logs directory not found"
fi

# Provide useful commands
echo -e "\n${YELLOW}Useful Commands:${NC}"
echo "  Cancel all array jobs:     scancel -n af3_msa_array -u $USER"
echo "  Cancel normal jobs only:   scancel -n af3_msa_array -u $USER -p normal"
echo "  Cancel hns jobs only:      scancel -n af3_msa_array -u $USER -p hns"
echo "  Cancel specific array:     scancel <array_job_id>"
echo "  View specific job output:  tail -f logs/<job_id>_<array_index>_MSA.out"
echo "  View job details:          scontrol show job <job_id>"
echo "  Clean up tmp files:        rm ${AF3_BASE_DIR}/msa_array_jobs_part*.tmp"
echo "  Resubmit failed jobs:      Check logs and rerun submit_msa_arrays.sh with updated CSV"

# Check for tmp files
echo -e "\n${BLUE}Temporary Files:${NC}"
if ls ${AF3_BASE_DIR}/msa_array_jobs_part*.tmp 1> /dev/null 2>&1; then
    echo "  Found partition files:"
    ls -lh ${AF3_BASE_DIR}/msa_array_jobs_part*.tmp | awk '{print "    " $9 " (" $5 ")"}'
else
    echo "  No partition files found"
fi

# Show progress estimation
RUNNING=$(count_jobs RUNNING)
PENDING=$(count_jobs PENDING)
TOTAL=$((RUNNING + PENDING))

if [ $TOTAL -gt 0 ]; then
    echo -e "\n${BLUE}Progress Estimation:${NC}"
    echo "  Jobs remaining: $TOTAL"
    if [ $RUNNING -gt 0 ]; then
        # Rough estimation assuming 2-4 hours per job
        MIN_HOURS=$(( (TOTAL * 2) / RUNNING ))
        MAX_HOURS=$(( (TOTAL * 4) / RUNNING ))
        echo "  Estimated time remaining: ${MIN_HOURS}-${MAX_HOURS} hours (rough estimate)"
    fi
fi
