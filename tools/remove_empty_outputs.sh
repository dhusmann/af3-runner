#!/bin/bash

# Base directory for jobs
BASE_DIR="/scratch/groups/ogozani/alphafold3/jobs/human_test_set"

# Counter variables
TOTAL=0
REMOVED=0
NOT_FOUND=0
NOT_EMPTY=0

# Read job names from stdin or from a file
while IFS= read -r job_name; do
    # Skip empty lines
    [ -z "$job_name" ] && continue
    
    TOTAL=$((TOTAL + 1))
    
    # Construct the output directory path
    output_dir="${BASE_DIR}/${job_name}/output"
    
    # Check if the output directory exists
    if [ -d "$output_dir" ]; then
        # Check if directory is empty (no files or subdirectories)
        if [ -z "$(ls -A "$output_dir")" ]; then
            # Directory is empty, remove it
            echo "Removing empty directory: $output_dir"
            rmdir "$output_dir"
            if [ $? -eq 0 ]; then
                REMOVED=$((REMOVED + 1))
            else
                echo "  ERROR: Failed to remove $output_dir"
            fi
        else
            echo "Skipping non-empty directory: $output_dir"
            NOT_EMPTY=$((NOT_EMPTY + 1))
        fi
    else
        echo "Directory not found: $output_dir"
        NOT_FOUND=$((NOT_FOUND + 1))
    fi
done

echo ""
echo "Summary:"
echo "  Total jobs processed: $TOTAL"
echo "  Empty directories removed: $REMOVED"
echo "  Non-empty directories skipped: $NOT_EMPTY"
echo "  Directories not found: $NOT_FOUND"
