#!/bin/bash

# Script to clean up temporary partition files after MSA jobs complete

AF3_BASE_DIR="/scratch/groups/ogozani/alphafold3"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}MSA Temporary File Cleanup${NC}"
echo "=========================="

# Check if any MSA array jobs are still running
RUNNING_JOBS=$(squeue -u $USER -n af3_msa_array -h | wc -l)

if [ $RUNNING_JOBS -gt 0 ]; then
    echo -e "${YELLOW}Warning: You still have $RUNNING_JOBS MSA array jobs running!${NC}"
    echo "Running jobs may still need these temporary files."
    echo
    squeue -u $USER -n af3_msa_array -o "%.18i %.9P %.50j %.8u %.2t %.10M %.6D %R" | head -10
    echo
    read -p "Are you sure you want to delete the temporary files? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cleanup cancelled."
        exit 0
    fi
fi

# Find and list temporary files
TMP_FILES=$(ls ${AF3_BASE_DIR}/msa_array_jobs_part*.tmp 2>/dev/null)

if [ -z "$TMP_FILES" ]; then
    echo "No temporary partition files found."
    exit 0
fi

echo "Found the following temporary files:"
for file in $TMP_FILES; do
    echo "  - $file ($(stat -c%s "$file" | numfmt --to=iec-i --suffix=B))"
done

echo
read -p "Delete these files? (y/N): " -n 1 -r
echo

if [[ $REPLY =~ ^[Yy]$ ]]; then
    for file in $TMP_FILES; do
        rm -v "$file"
    done
    echo -e "${GREEN}Cleanup complete!${NC}"
else
    echo "Cleanup cancelled."
fi
