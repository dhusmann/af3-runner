#!/bin/bash

# sync_all.sh - Complete AlphaFold3 sync workflow orchestrator
# 
# This script runs the complete parallel sync workflow:
# 1. Submits parallel rsync, MSA archiving, and seed compression jobs
# 2. Submits Google Drive sync job for outputs, scripts, and analysis
# 
# All operations run on compute nodes via SLURM - no login node blocking.

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Configuration
DRY_RUN=0
SKIP_GDRIVE=0
QUIET=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--dry-run)
            DRY_RUN=1
            shift
            ;;
        --local-only)
            SKIP_GDRIVE=1
            shift
            ;;
        -q|--quiet)
            QUIET=1
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo
            echo "Complete AlphaFold3 sync workflow - submits all sync jobs to SLURM"
            echo
            echo "Options:"
            echo "  -n, --dry-run     Preview jobs without submitting"
            echo "  -q, --quiet       Minimal output"
            echo "  --local-only      Skip Google Drive sync (local processing only)"
            echo "  -h, --help        Show this help message"
            echo
            echo "This script submits:"
            echo "  • Parallel rsync job (organizes outputs on compute nodes)"
            echo "  • MSA archiving job (creates deduplicated archives)"
            echo "  • Seed compression job (compresses seed directories)"
            echo "  • Google Drive sync job (uploads everything including scripts)"
            echo
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h for help"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}AlphaFold3 Complete Sync Workflow${NC}"
echo "=================================="
echo "Automated parallel processing on compute nodes"
echo

# Build options for subscripts
SYNC_OPTS=""
if [ $DRY_RUN -eq 1 ]; then
    SYNC_OPTS="--dry-run"
    echo -e "${YELLOW}DRY RUN MODE - No jobs will be submitted${NC}"
    echo
fi
if [ $QUIET -eq 1 ]; then
    SYNC_OPTS="$SYNC_OPTS --quiet"
fi

# Step 1: Submit local sync jobs (rsync, MSA archiving, seed compression)
echo -e "${GREEN}Step 1: Submitting local sync and processing jobs${NC}"
echo "This will submit parallel jobs for:"
echo "  • Rsync operations (organize outputs)"
echo "  • MSA archiving (deduplicated archives)"  
echo "  • Seed compression (compress seed directories)"
echo

./sync_organize_outputs.sh $SYNC_OPTS

if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to submit local sync jobs${NC}"
    exit 1
fi

if [ $DRY_RUN -eq 1 ]; then
    echo
    echo -e "${YELLOW}Dry run mode - no actual jobs submitted${NC}"
    if [ $SKIP_GDRIVE -eq 0 ]; then
        echo "Would also submit Google Drive sync job in real run"
    fi
    echo "Run without -n flag to submit all jobs"
    exit 0
fi

echo
echo -e "${GREEN}Local sync jobs submitted successfully${NC}"

# Step 2: Submit Google Drive sync job (unless skipped)
if [ $SKIP_GDRIVE -eq 1 ]; then
    echo -e "${YELLOW}Skipping Google Drive sync (--local-only specified)${NC}"
else
    echo
    echo -e "${GREEN}Step 2: Submitting Google Drive sync job${NC}"
    echo "This will backup and sync:"
    echo "  • All scripts and configuration files"
    echo "  • Tools directory"
    echo "  • Analysis directory"
    echo "  • All organized output data"
    echo

    ./rclone_to_gdrive.sh

    if [ $? -ne 0 ]; then
        echo -e "${RED}Failed to submit Google Drive sync job${NC}"
        echo "Local sync jobs are still running - check with: squeue -u \$USER"
        exit 1
    fi
    
    echo
    echo -e "${GREEN}Google Drive sync job submitted successfully${NC}"
fi

# Display monitoring information
echo
echo -e "${GREEN}All sync jobs submitted!${NC}"
echo "========================="

# Show active jobs
if ! $QUIET; then
    echo
    echo -e "${BLUE}Active jobs:${NC}"
    if command -v squeue >/dev/null 2>&1; then
        squeue -u $USER -o "%.18i %.12j %.8T %.10M %.6D %R" 2>/dev/null || echo "  Run: squeue -u \$USER"
    else
        echo "  Run: squeue -u \$USER"
    fi
fi

echo
echo -e "${BLUE}Monitoring commands:${NC}"
echo "==================="
echo "Check all jobs:      squeue -u \$USER"
echo "Check job details:   squeue -j JOB_ID"
echo "View logs:"
echo "  Rsync logs:        ls logs/sync_rsync_*"
echo "  MSA archiving:     ls logs/archive_msa_*"  
echo "  Seed compression:  ls logs/pack-seeds_*"
if [ $SKIP_GDRIVE -eq 0 ]; then
    echo "  Google Drive sync: ls logs/Rclone_*"
fi
echo
echo "Live monitoring:"
echo "  tail -f logs/sync_rsync_*.out      # Rsync progress"
echo "  tail -f logs/archive_msa_*.out     # MSA archiving"
echo "  tail -f logs/pack-seeds_*.out      # Seed compression"
if [ $SKIP_GDRIVE -eq 0 ]; then
    echo "  tail -f logs/Rclone_*.out          # Google Drive sync"
fi

echo
echo -e "${GREEN}Workflow complete!${NC} All operations are running on compute nodes."
echo "Your login node is free for other work."

if [ $SKIP_GDRIVE -eq 0 ]; then
    echo
    echo "📁 After completion, your data will be organized in Google Drive at:"
    echo "   • gozani_labshare_alphafold:alphafold3/scripts/     (your scripts)"
    echo "   • gozani_labshare_alphafold:alphafold3/analysis/    (analysis files)" 
    echo "   • gozani_labshare_alphafold:alphafold3/output/      (organized results)"
fi