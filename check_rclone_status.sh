#!/bin/bash

# check_rclone_status.sh - Check status of rclone sync jobs
# Shows active, completed, and failed rclone jobs

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}Rclone Sync Job Status${NC}"
echo "======================"
echo

# Check active rclone jobs
echo -e "${GREEN}Active sync jobs:${NC}"
ACTIVE_JOBS=$(squeue -u $USER -n "SyncAF3_Gdrive,SyncAF3_Retry" -h)

if [ -z "$ACTIVE_JOBS" ]; then
    echo "  No active rclone jobs"
else
    echo "$ACTIVE_JOBS" | while read -r line; do
        JOB_ID=$(echo "$line" | awk '{print $1}')
        JOB_NAME=$(echo "$line" | awk '{print $3}')
        STATE=$(echo "$line" | awk '{print $5}')
        TIME=$(echo "$line" | awk '{print $6}')
        
        echo "  Job $JOB_ID ($JOB_NAME): $STATE for $TIME"
    done
fi

echo

# Check recent completed jobs
echo -e "${GREEN}Recent rclone jobs (last 24h):${NC}"
RECENT_JOBS=$(sacct -u $USER -S $(date -d '24 hours ago' +%Y-%m-%d) -n \
    --name="SyncAF3_Gdrive,SyncAF3_Retry" \
    --format=JobID,JobName%20,State,ExitCode,Start,Elapsed)

if [ -z "$RECENT_JOBS" ]; then
    echo "  No recent jobs"
else
    echo "$RECENT_JOBS" | grep -v ".batch" | while read -r line; do
        JOB_ID=$(echo "$line" | awk '{print $1}')
        STATE=$(echo "$line" | awk '{print $3}')
        EXIT_CODE=$(echo "$line" | awk '{print $4}')
        
        if [[ "$STATE" == "COMPLETED" ]]; then
            echo -e "  ${GREEN}✓${NC} Job $JOB_ID: $STATE"
        elif [[ "$STATE" == "FAILED" ]] || [[ "$EXIT_CODE" != "0:0" ]]; then
            echo -e "  ${RED}✗${NC} Job $JOB_ID: $STATE (exit: $EXIT_CODE)"
        else
            echo -e "  ${YELLOW}?${NC} Job $JOB_ID: $STATE"
        fi
    done
fi

echo

# Check latest log file
LATEST_LOG=$(ls -t logs/Rclone_*.out 2>/dev/null | head -1)

if [ -f "$LATEST_LOG" ]; then
    echo -e "${GREEN}Latest log file:${NC} $LATEST_LOG"
    
    # Check for quota errors in recent log
    if grep -q "quotaExceeded\|userRateLimitExceeded" "$LATEST_LOG" 2>/dev/null; then
        echo -e "${YELLOW}  Warning: Quota limit detected in recent sync${NC}"
        LAST_QUOTA=$(grep -E "quotaExceeded|userRateLimitExceeded" "$LATEST_LOG" | tail -1)
        echo "  Last quota error: $(echo "$LAST_QUOTA" | cut -c1-60)..."
    fi
    
    # Check for successful transfers
    TRANSFERRED=$(grep -c "Copied (new)" "$LATEST_LOG" 2>/dev/null || echo "0")
    if [ "$TRANSFERRED" -gt 0 ]; then
        echo "  Files transferred: $TRANSFERRED"
    fi
    
    # Show last few lines
    echo
    echo "  Last 5 lines:"
    tail -5 "$LATEST_LOG" | sed 's/^/    /'
else
    echo "No rclone log files found in logs/"
fi

echo

# Check if retry is scheduled
PENDING_RETRY=$(squeue -u $USER -n "SyncAF3_Retry" -t PENDING -h)
if [ -n "$PENDING_RETRY" ]; then
    echo -e "${YELLOW}Scheduled retry job:${NC}"
    echo "$PENDING_RETRY" | while read -r line; do
        JOB_ID=$(echo "$line" | awk '{print $1}')
        START_TIME=$(squeue -j $JOB_ID -o "%S" -h)
        echo "  Job $JOB_ID scheduled for: $START_TIME"
    done
    echo
fi

# Quick stats
echo -e "${BLUE}Quick Actions:${NC}"
echo "  View active job output:  tail -f logs/Rclone_<job_id>.out"
echo "  Cancel sync job:         scancel <job_id>"
echo "  Submit new sync:         ./rclone_to_gdrive.sh"
echo "  Check what's ready:      ./check_sync_status.sh"
