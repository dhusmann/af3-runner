#!/bin/bash

# pipeline_status.sh - Quick dashboard view of the AlphaFold3 pipeline
# Shows current state of all components at a glance

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${BOLD}AlphaFold3 Pipeline Status Dashboard${NC}"
echo "===================================="
date
echo

# Check if pipeline is running
CYCLE_JOBS=$(squeue -u $USER -n af3_cycle -h | wc -l)
if [ $CYCLE_JOBS -gt 0 ]; then
    CYCLE_INFO=$(squeue -u $USER -n af3_cycle -h -o "%i %S %j" | head -1)
    echo -e "${GREEN}● Pipeline Status: RUNNING${NC}"
    echo "  $CYCLE_INFO"
    
    # Try to determine current cycle from job name or output
    LATEST_LOG=$(ls -t af3_cycle_*.out 2>/dev/null | head -1)
    if [ -f "$LATEST_LOG" ]; then
        CURRENT_CYCLE=$(grep "AlphaFold3 Cycle" "$LATEST_LOG" | tail -1 | grep -oP 'Cycle \K\d+')
        if [ -n "$CURRENT_CYCLE" ]; then
            echo "  Current cycle: $CURRENT_CYCLE of 24"
        fi
    fi
else
    echo -e "${YELLOW}● Pipeline Status: NOT RUNNING${NC}"
fi

echo

# MSA Array Jobs
echo -e "${BOLD}MSA Generation:${NC}"
MSA_RUNNING=$(squeue -u $USER -n af3_msa_array -t RUNNING -h | wc -l)
MSA_PENDING=$(squeue -u $USER -n af3_msa_array -t PENDING -h | wc -l)
MSA_TOTAL=$((MSA_RUNNING + MSA_PENDING))

if [ $MSA_TOTAL -gt 0 ]; then
    echo -e "  ${BLUE}● MSA Arrays: $MSA_RUNNING running, $MSA_PENDING pending${NC}"
    
    # Show partition breakdown
    MSA_NORMAL_R=$(squeue -u $USER -n af3_msa_array -p normal -t RUNNING -h | wc -l)
    MSA_HNS_R=$(squeue -u $USER -n af3_msa_array -p hns -t RUNNING -h | wc -l)
    echo "    Normal partition: $MSA_NORMAL_R running"
    echo "    HNS partition: $MSA_HNS_R running"
else
    echo -e "  ${GREEN}✓ No MSA jobs running${NC}"
fi

# Check MSA completion rate if file exists
if [ -f "msa_array_jobs.csv" ]; then
    TOTAL_MSA=$(tail -n +2 msa_array_jobs.csv | grep -v "^$" | wc -l)
    echo "  Total MSA jobs in queue: $TOTAL_MSA"
fi

echo

# GPU Jobs
echo -e "${BOLD}GPU Inference:${NC}"
GPU_RUNNING=$(squeue -u $USER -p gpu -h | wc -l)
echo -e "  ${CYAN}● GPU Jobs: $GPU_RUNNING running${NC}"

echo

# File Status
echo -e "${BOLD}Input Files:${NC}"
for file in folding_jobs.csv msa_array_jobs.csv waiting_for_msa.csv; do
    if [ -f "$file" ]; then
        COUNT=$(tail -n +2 "$file" | grep -v "^$" | wc -l)
        echo -e "  ${GREEN}✓${NC} $file ($COUNT jobs)"
    else
        echo -e "  ${RED}✗${NC} $file"
    fi
done

echo

# Quick job stage summary
echo -e "${BOLD}Job Progress Summary:${NC}"
if [ -f "folding_jobs.csv" ]; then
    # Quick counts without full scan
    STAGE1=0
    STAGE2=0
    STAGE3=0
    
    # Sample first 100 jobs for quick estimate
    SAMPLE_SIZE=100
    JOBS=$(tail -n +2 folding_jobs.csv | grep -v "^$" | head -$SAMPLE_SIZE | awk '{print $1}')
    TOTAL_SAMPLED=0
    
    while IFS= read -r job_name; do
        if [ -n "$job_name" ]; then
            ((TOTAL_SAMPLED++))
            
            # Quick check without detailed scan
            if [ -d "jobs/$job_name/output" ] || [ -d "jobs/human_test_set/$job_name/output" ]; then
                ((STAGE3++))
            elif [ -d "jobs/$job_name/output_msa" ] || [ -d "jobs/human_test_set/$job_name/output_msa" ]; then
                if ls jobs/$job_name/output_msa/*.json >/dev/null 2>&1 || ls jobs/human_test_set/$job_name/output_msa/*.json >/dev/null 2>&1; then
                    ((STAGE2++))
                else
                    ((STAGE1++))
                fi
            else
                ((STAGE1++))
            fi
        fi
    done <<< "$JOBS"
    
    if [ $TOTAL_SAMPLED -gt 0 ]; then
        echo "  Based on sample of $TOTAL_SAMPLED jobs:"
        echo -e "  ${YELLOW}Stage 1 (Need MSA):${NC} ~$((STAGE1 * 100 / TOTAL_SAMPLED))%"
        echo -e "  ${BLUE}Stage 2 (Need GPU):${NC} ~$((STAGE2 * 100 / TOTAL_SAMPLED))%"
        echo -e "  ${GREEN}Stage 3 (Complete):${NC} ~$((STAGE3 * 100 / TOTAL_SAMPLED))%"
    fi
fi

echo

# Recent completions
echo -e "${BOLD}Recent Activity:${NC}"
if [ -f "logged_folding_jobs.csv" ]; then
    RECENT=$(tail -5 logged_folding_jobs.csv | grep -E "(GPU|MSA).*,0$" | tail -3)
    if [ -n "$RECENT" ]; then
        echo "  Recent completions:"
        echo "$RECENT" | while IFS=',' read -r job_name job_id stage start end duration status; do
            echo "    $job_name - $stage completed in $duration"
        done
    else
        echo "  No recent completions"
    fi
else
    echo "  No log file found"
fi

echo
echo -e "${CYAN}For detailed status: ./get_job_status.sh${NC}"
echo -e "${CYAN}For MSA monitoring:  ./monitor_msa_arrays.sh${NC}"
