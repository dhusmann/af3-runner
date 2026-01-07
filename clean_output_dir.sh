#!/bin/bash

# clean_output_dir.sh - Clean up the output directory
# Useful for starting fresh or removing old syncs

# Colors
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'

OUTPUT_DIR="/scratch/groups/ogozani/alphafold3/output"

echo -e "${YELLOW}AlphaFold3 Output Directory Cleanup${NC}"
echo "==================================="
echo

# Check if output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory does not exist: $OUTPUT_DIR"
    echo "Nothing to clean."
    exit 0
fi

# Show current status
echo "Current output directory contents:"
echo

# Count files and show size
TOTAL_FILES=$(find "$OUTPUT_DIR" -type f | wc -l)
TOTAL_SIZE=$(du -sh "$OUTPUT_DIR" 2>/dev/null | cut -f1)

echo "Total files: $TOTAL_FILES"
echo "Total size: $TOTAL_SIZE"
echo

# Show directory structure
echo "Directory structure:"
for dir in "$OUTPUT_DIR"/*/; do
    if [ -d "$dir" ]; then
        dir_name=$(basename "$dir")
        file_count=$(find "$dir" -type f | wc -l)
        dir_size=$(du -sh "$dir" 2>/dev/null | cut -f1)
        echo "  $dir_name: $file_count files ($dir_size)"
    fi
done

echo
echo -e "${RED}WARNING: This will delete all files in the output directory!${NC}"
echo "This action cannot be undone."
echo

read -p "Are you sure you want to delete $OUTPUT_DIR? Type 'yes' to confirm: " -r
echo

if [ "$REPLY" = "yes" ]; then
    echo "Removing output directory..."
    rm -rf "$OUTPUT_DIR"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Output directory removed successfully${NC}"
        echo
        echo "You can now run sync_organize_outputs.sh to create a fresh sync"
    else
        echo -e "${RED}Failed to remove output directory${NC}"
        echo "Check permissions and try again"
        exit 1
    fi
else
    echo "Cleanup cancelled"
    echo
    echo "To remove specific subdirectories instead:"
    echo "  rm -rf $OUTPUT_DIR/project_name"
fi
