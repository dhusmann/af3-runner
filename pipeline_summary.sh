#!/bin/bash

# pipeline_summary.sh - Complete pipeline summary including sync status
# Provides a comprehensive overview of the entire pipeline state

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# Configuration
BASE_DIR="/scratch/groups/ogozani/alphafold3"

echo -e "${BOLD}${BLUE}AlphaFold3 Pipeline Complete Summary${NC}"
echo "===================================="
date
echo

# Section 1: Pipeline Status
echo -e "${BOLD}1. Pipeline Automation${NC}"
CYCLE_JOBS=$(squeue -u $USER -n af3_cycle -h | wc -l)
if [ $CYCLE_JOBS -gt 0 ]; then
    echo -e "   ${GREEN}● Status: RUNNING${NC}"
    LATEST_LOG=$(ls -t af3_cycle_*.out 2>/dev/null | head -1)
    if [ -f "$LATEST_LOG" ]; then
        CURRENT_CYCLE=$(grep "AlphaFold3 Cycle" "$LATEST_LOG" | tail -1 | grep -oP 'Cycle \K\d+')
        if [ -n "$CURRENT_CYCLE" ]; then
            echo "   Current cycle: $CURRENT_CYCLE of 24"
            PERCENT=$((CURRENT_CYCLE * 100 / 24))
            echo "   Progress: ${PERCENT}%"
        fi
    fi
else
    echo -e "   ${YELLOW}● Status: NOT RUNNING${NC}"
fi
echo

# Section 2: Active Jobs
echo -e "${BOLD}2. Active Jobs${NC}"
MSA_RUNNING=$(squeue -u $USER -n af3_msa_array -t RUNNING -h | wc -l)
MSA_PENDING=$(squeue -u $USER -n af3_msa_array -t PENDING -h | wc -l)
GPU_RUNNING=$(squeue -u $USER -p gpu -h | wc -l)

echo "   MSA Arrays: $MSA_RUNNING running, $MSA_PENDING pending"
echo "   GPU Jobs: $GPU_RUNNING running"
echo

# Section 3: Job Progress
echo -e "${BOLD}3. Job Progress${NC}"
if [ -f "folding_jobs.csv" ]; then
    TOTAL_JOBS=$(tail -n +2 folding_jobs.csv | grep -v "^$" | wc -l)
    
    # Quick sample check
    SAMPLE_SIZE=50
    JOBS=$(tail -n +2 folding_jobs.csv | grep -v "^$" | head -$SAMPLE_SIZE | awk '{print $1}')
    STAGE1=0
    STAGE2=0
    STAGE3=0
    
    while IFS= read -r job_name; do
        if [ -n "$job_name" ]; then
            if [ -d "$BASE_DIR/jobs/$job_name/output" ] || [ -d "$BASE_DIR/jobs/human_test_set/$job_name/output" ]; then
                if ls "$BASE_DIR/jobs/$job_name/output/"*/*.cif >/dev/null 2>&1 || ls "$BASE_DIR/jobs/human_test_set/$job_name/output/"*/*.cif >/dev/null 2>&1; then
                    ((STAGE3++))
                fi
            elif [ -d "$BASE_DIR/jobs/$job_name/output_msa" ] || [ -d "$BASE_DIR/jobs/human_test_set/$job_name/output_msa" ]; then
                ((STAGE2++))
            else
                ((STAGE1++))
            fi
        fi
    done <<< "$JOBS"
    
    # Extrapolate
    if [ $SAMPLE_SIZE -lt $TOTAL_JOBS ]; then
        STAGE1_EST=$((STAGE1 * TOTAL_JOBS / SAMPLE_SIZE))
        STAGE2_EST=$((STAGE2 * TOTAL_JOBS / SAMPLE_SIZE))
        STAGE3_EST=$((STAGE3 * TOTAL_JOBS / SAMPLE_SIZE))
        echo "   Total jobs: $TOTAL_JOBS"
        echo "   Estimated breakdown (based on sample):"
        echo -e "   ${YELLOW}Need MSA:${NC} ~$STAGE1_EST jobs"
        echo -e "   ${BLUE}Need GPU:${NC} ~$STAGE2_EST jobs"
        echo -e "   ${GREEN}Complete:${NC} ~$STAGE3_EST jobs"
        
        if [ $TOTAL_JOBS -gt 0 ]; then
            PERCENT_COMPLETE=$((STAGE3_EST * 100 / TOTAL_JOBS))
            echo "   Overall: ~${PERCENT_COMPLETE}% complete"
        fi
    else
        echo "   Total jobs: $TOTAL_JOBS"
        echo -e "   ${YELLOW}Need MSA:${NC} $STAGE1 jobs"
        echo -e "   ${BLUE}Need GPU:${NC} $STAGE2 jobs"
        echo -e "   ${GREEN}Complete:${NC} $STAGE3 jobs"
    fi
else
    echo "   No folding_jobs.csv found"
fi
echo

# Section 4: Output Sync Status
echo -e "${BOLD}4. Output Sync Status${NC}"

# Count ready outputs
READY_COUNT=0
for job_dir in "$BASE_DIR/jobs/"*/ "$BASE_DIR/jobs/human_test_set/"*/; do
    if [ -d "$job_dir" ] && [ -d "${job_dir}output" ] && [ -n "$(ls -A "${job_dir}output" 2>/dev/null)" ]; then
        ((READY_COUNT++))
    fi
done

echo "   Outputs ready to sync: $READY_COUNT"

# Check if already synced
if [ -d "$BASE_DIR/output" ]; then
    SYNCED_COUNT=$(find "$BASE_DIR/output" -maxdepth 2 -type d -name "*" | grep -v "^$BASE_DIR/output$" | wc -l)
    echo "   Already synced: $SYNCED_COUNT directories"
fi

# Estimate size
if [ $READY_COUNT -gt 0 ]; then
    SYNC_SIZE=$(find "$BASE_DIR/jobs" -path "*/output/*" -type f ! -name "*_data.json" ! -name "TERMS_OF_USE.md" -exec du -ch {} + 2>/dev/null | grep total$ | head -1 | cut -f1)
    if [ -n "$SYNC_SIZE" ]; then
        echo "   Estimated sync size: $SYNC_SIZE (excluding *_data.json)"
    fi
fi

# Check rclone
if command -v rclone &> /dev/null && rclone listremotes | grep -q "^gdrive:$"; then
    echo -e "   ${GREEN}✓ Rclone configured${NC}"
else
    echo -e "   ${YELLOW}✗ Rclone not configured${NC} (run ./setup_rclone_gdrive.sh)"
fi
echo

# Section 5: Recent Activity
echo -e "${BOLD}5. Recent Activity${NC}"
if [ -f "logged_folding_jobs.csv" ]; then
    RECENT=$(tail -10 logged_folding_jobs.csv | grep -E "(GPU|MSA).*,0$" | tail -5)
    if [ -n "$RECENT" ]; then
        echo "$RECENT" | while IFS=',' read -r job_name job_id stage start end duration status; do
            echo "   $job_name - $stage completed in $duration"
        done
    else
        echo "   No recent completions"
    fi
else
    echo "   No activity log found"
fi
echo

# Section 6: Next Steps
echo -e "${BOLD}6. Recommended Next Steps${NC}"

if [ $CYCLE_JOBS -eq 0 ] && [ $READY_COUNT -eq 0 ]; then
    echo "   → Start pipeline: ./launch_af3.sh"
elif [ $READY_COUNT -gt 10 ]; then
    echo "   → Sync outputs: ./sync_all.sh"
    echo "   → Check details: ./get_job_status.sh"
elif [ $GPU_RUNNING -eq 0 ] && [ $STAGE2 -gt 0 ]; then
    echo "   → Submit GPU jobs: ./submit_dist.sh"
fi

if [ $MSA_RUNNING -eq 0 ] && [ -f "msa_array_jobs.csv" ]; then
    MSA_TODO=$(tail -n +2 msa_array_jobs.csv | grep -v "^$" | wc -l)
    if [ $MSA_TODO -gt 0 ]; then
        echo "   → Submit MSA jobs: ./submit_msa_arrays.sh"
    fi
fi

echo
echo -e "${CYAN}For detailed information, use:${NC}"
echo "   ./get_job_status.sh -v      # Detailed job status"
echo "   ./check_sync_status.sh      # Sync readiness"
echo "   ./monitor_msa_arrays.sh     # MSA job details"
