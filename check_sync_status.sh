#!/bin/bash

# check_sync_status.sh - Check which outputs are ready to sync
# Shows what would be synced without actually doing it

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

BASE_DIR="/scratch/groups/ogozani/alphafold3"
JOBS_DIR="${BASE_DIR}/jobs"

echo -e "${BLUE}AlphaFold3 Output Sync Status${NC}"
echo "=============================="
echo

# Count jobs with outputs
TOTAL_WITH_OUTPUT=0
TOTAL_SIZE=0

# Check top-level jobs
echo "Checking jobs with completed outputs..."
echo

echo "Top-level jobs:"
for job_dir in "${JOBS_DIR}"/*/; do
    if [ -d "$job_dir" ] && [ "$(basename "$job_dir")" != "human_test_set" ]; then
        job_name=$(basename "$job_dir")
        output_dir="${job_dir}output/"
        
        if [ -d "$output_dir" ] && [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then
            # Get size excluding *_data.json
            size=$(find "$output_dir" -type f ! -name "*_data.json" ! -name "TERMS_OF_USE.md" -exec du -ch {} + | grep total$ | cut -f1)
            echo "  ✓ $job_name ($size)"
            ((TOTAL_WITH_OUTPUT++))
        fi
    fi
done

echo
echo "Human test set jobs:"
if [ -d "${JOBS_DIR}/human_test_set" ]; then
    for job_dir in "${JOBS_DIR}/human_test_set"/*/; do
        if [ -d "$job_dir" ]; then
            job_name=$(basename "$job_dir")
            output_dir="${job_dir}output/"
            
            if [ -d "$output_dir" ] && [ -n "$(ls -A "$output_dir" 2>/dev/null)" ]; then
                # Get size excluding *_data.json
                size=$(find "$output_dir" -type f ! -name "*_data.json" ! -name "TERMS_OF_USE.md" -exec du -ch {} + | grep total$ | cut -f1)
                echo "  ✓ $job_name ($size)"
                ((TOTAL_WITH_OUTPUT++))
            fi
        fi
    done
fi

echo
echo "Summary:"
echo "--------"
echo "Jobs with outputs ready to sync: $TOTAL_WITH_OUTPUT"

# Calculate total size to sync
echo
echo "Calculating total size to sync (excluding *_data.json)..."
SYNC_SIZE=$(find "${JOBS_DIR}" -path "*/output/*" -type f ! -name "*_data.json" ! -name "TERMS_OF_USE.md" -exec du -ch {} + 2>/dev/null | grep total$ | cut -f1)
EXCLUDED_SIZE=$(find "${JOBS_DIR}" -path "*/output/*" -name "*_data.json" -exec du -ch {} + 2>/dev/null | grep total$ | cut -f1)

echo -e "${GREEN}Size to sync: $SYNC_SIZE${NC}"
if [ -n "$EXCLUDED_SIZE" ]; then
    echo -e "${YELLOW}Size excluded (*_data.json): $EXCLUDED_SIZE${NC}"
fi

# Check if output directory exists
echo
if [ -d "${BASE_DIR}/output" ]; then
    ALREADY_SYNCED=$(find "${BASE_DIR}/output" -type f | wc -l)
    if [ $ALREADY_SYNCED -gt 0 ]; then
        echo -e "${BLUE}Note: ${BASE_DIR}/output already contains $ALREADY_SYNCED files${NC}"
        echo "Running sync_organize_outputs.sh will update with any new/changed files"
    fi
else
    echo "Output directory does not exist yet"
    echo "Run sync_organize_outputs.sh to create and populate it"
fi

echo
echo "Next step: ./sync_organize_outputs.sh"
