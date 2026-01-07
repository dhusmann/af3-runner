#!/bin/bash

# archive_msa_data.sh - Archive unique AlphaFold MSA data files
# 
# This script creates compressed archives of unique MSA data files,
# avoiding duplicate MSAs that were reused across multiple jobs.
#
# Usage:
#   ./archive_msa_data.sh [job_group]
#
# If job_group is specified, only that group will be processed.
# If no job_group is specified, all unarchived groups will be processed.
#
# Environment variables:
#   COMPRESSION_THREADS (default 4) - threads for pigz compression
#   DRY_RUN (default 0) - set to 1 to preview without creating archives

set -euo pipefail

# Configuration
BASE_DIR="/scratch/groups/ogozani/alphafold3"
JOBS_DIR="${BASE_DIR}/jobs"
OUTPUT_MSA_DIR="${BASE_DIR}/output/msa"
MASTER_INDEX="${OUTPUT_MSA_DIR}/master_index.csv"
COMPRESSION_THREADS="${COMPRESSION_THREADS:-4}"
DRY_RUN="${DRY_RUN:-0}"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Statistics
TOTAL_UNIQUE_FILES=0
TOTAL_SKIPPED_FILES=0
TOTAL_SPACE_SAVED=0
TOTAL_ARCHIVES_CREATED=0

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_MSA_DIR"

# Initialize master index if it doesn't exist
if [[ ! -f "$MASTER_INDEX" ]]; then
    echo "archive_name,job_directory,json_filename,file_size_bytes,date_archived" > "$MASTER_INDEX"
fi

# Function to extract job group from directory name
get_job_group() {
    local job_dir="$1"
    # Extract everything before the first hyphen
    echo "${job_dir%%-*}"
}

# Function to check if file is already archived
is_already_archived() {
    local job_dir="$1"
    local json_file="$2"
    
    # Simple grep check against master index
    if grep -q ",$job_dir,$json_file," "$MASTER_INDEX" 2>/dev/null; then
        return 0
    fi
    
    return 1
}

# Function to find unique MSA files for a job group
find_unique_msa_files() {
    local job_group="$1"
    local temp_list=$(mktemp)
    
    echo -e "${BLUE}Finding unique MSA files for job group: $job_group${NC}" >&2
    
    # Find job directories first, then look for MSA files in each
    # This is more efficient than using path filtering with find
    for jobs_root in "$JOBS_DIR" "$JOBS_DIR/human_test_set"; do
        if [[ -d "$jobs_root" ]]; then
            for job_path in "$jobs_root"/${job_group}-*; do
                if [[ -d "$job_path" ]]; then
                    job_name=$(basename "$job_path")
                    msa_dir="${job_path}/output_msa"
                    
                    if [[ -d "$msa_dir" ]]; then
                        # Look for unique MSA files in this directory
                        for json_file in "$msa_dir"/*_data.json "$msa_dir"/*/*_data.json; do
                            if [[ -f "$json_file" && "$(basename "$json_file")" != "alphafold_input_with_msa.json" ]]; then
                                json_basename=$(basename "$json_file")
                                
                                # Check if already archived
                                if ! is_already_archived "$job_name" "$json_basename"; then
                                    file_size=$(stat -c%s "$json_file" 2>/dev/null || echo 0)
                                    echo "${job_path}|${json_file}|${file_size}" >> "$temp_list"
                                fi
                            fi
                        done
                    fi
                fi
            done
        fi
    done
    echo "$temp_list"
}

# Function to create archive for a job group
archive_job_group() {
    local job_group="$1"
    local file_list="$2"
    
    # Check if there are files to archive
    local file_count=$(wc -l < "$file_list" 2>/dev/null || echo 0)
    if [[ $file_count -eq 0 ]]; then
        echo "  No unique MSA files found for group: $job_group (all likely use shared MSAs)"
        return 0
    fi
    
    # Generate archive name with date
    local date_str=$(date '+%Y_%m_%d')
    local archive_base="${job_group}_${date_str}"
    local archive_name="${archive_base}.tar.gz"
    local archive_path="${OUTPUT_MSA_DIR}/${archive_name}"
    
    # Handle duplicate archive names
    local counter=1
    while [[ -f "$archive_path" ]]; do
        archive_name="${archive_base}_${counter}.tar.gz"
        archive_path="${OUTPUT_MSA_DIR}/${archive_name}"
        ((counter++))
    done
    
    echo -e "${GREEN}Creating archive: $archive_name${NC}"
    echo "  Files to archive: $file_count"
    
    if [[ $DRY_RUN -eq 1 ]]; then
        echo -e "${YELLOW}  [DRY RUN] Would create $archive_path${NC}"
        while IFS='|' read -r job_path json_file file_size; do
            echo "    Would include: $(basename "$job_path")/${json_file##*/}"
        done < "$file_list"
        return 0
    fi
    
    # Create temporary staging directory
    local staging_dir=$(mktemp -d)
    local total_size=0
    
    echo "  Copying files to staging directory..."
    while IFS='|' read -r job_path json_file file_size; do
        job_name=$(basename "$job_path")
        json_basename=$(basename "$json_file")
        
        # Create job subdirectory in staging
        mkdir -p "${staging_dir}/${job_name}"
        
        # Copy file to staging
        cp "$json_file" "${staging_dir}/${job_name}/${json_basename}"
        
        total_size=$((total_size + file_size))
        echo "    Staged: ${job_name}/${json_basename} ($(numfmt --to=iec "$file_size"))"
    done < "$file_list"
    
    echo "  Total size: $(numfmt --to=iec "$total_size")"
    echo "  Creating compressed archive..."
    
    # Create compressed archive
    pushd "$staging_dir" >/dev/null
    if command -v pigz >/dev/null 2>&1; then
        tar -cf - */ | pigz -p "$COMPRESSION_THREADS" -7 > "$archive_path"
    else
        tar -czf "$archive_path" */
    fi
    popd >/dev/null
    
    # Verify archive
    if tar -tzf "$archive_path" >/dev/null 2>&1; then
        local archive_size=$(stat -c%s "$archive_path")
        local compression_ratio=$(echo "scale=1; $archive_size * 100 / $total_size" | bc 2>/dev/null || echo "N/A")
        echo -e "${GREEN}  ✓ Archive created successfully${NC}"
        echo "    Archive size: $(numfmt --to=iec "$archive_size")"
        echo "    Compression ratio: ${compression_ratio}%"
        
        # Update master index
        while IFS='|' read -r job_path json_file file_size; do
            job_name=$(basename "$job_path")
            json_basename=$(basename "$json_file")
            echo "${archive_name},${job_name},${json_basename},${file_size},$(date '+%Y-%m-%d')" >> "$MASTER_INDEX"
        done < "$file_list"
        
        # Update statistics
        TOTAL_UNIQUE_FILES=$((TOTAL_UNIQUE_FILES + file_count))
        TOTAL_ARCHIVES_CREATED=$((TOTAL_ARCHIVES_CREATED + 1))
        TOTAL_SPACE_SAVED=$((TOTAL_SPACE_SAVED + total_size - archive_size))
    else
        echo -e "${RED}  ✗ Archive verification failed${NC}"
        rm -f "$archive_path"
        cleanup_staging "$staging_dir"
        return 1
    fi
    
    # Cleanup staging directory
    cleanup_staging "$staging_dir"
}

# Function to cleanup staging directory
cleanup_staging() {
    local staging_dir="$1"
    if [[ -d "$staging_dir" ]]; then
        rm -rf "$staging_dir"
    fi
}

# Function to get all unique job groups
get_all_job_groups() {
    local groups_file=$(mktemp)
    
    # Find all job directories in both main and human_test_set
    {
        # Main jobs directory
        find "$JOBS_DIR" -maxdepth 1 -type d -name "*-*" | while read job_path; do
            job_name=$(basename "$job_path")
            if [[ "$job_name" != "human_test_set" ]]; then
                get_job_group "$job_name"
            fi
        done
        
        # Human test set subdirectory
        if [[ -d "$JOBS_DIR/human_test_set" ]]; then
            find "$JOBS_DIR/human_test_set" -maxdepth 1 -type d -name "*-*" | while read job_path; do
                job_name=$(basename "$job_path")
                get_job_group "$job_name"
            done
        fi
    } | sort -u > "$groups_file"
    
    echo "$groups_file"
}

# Function to count shared MSAs for statistics
count_shared_msas() {
    local job_group="$1"
    
    # Use find with -exec to count shared MSAs directly
    find "$JOBS_DIR" "$JOBS_DIR/human_test_set" -maxdepth 2 -name "${job_group}-*" -type d 2>/dev/null | \
    while IFS= read -r job_path; do
        msa_dir="${job_path}/output_msa"
        if [[ -d "$msa_dir" && -f "${msa_dir}/alphafold_input_with_msa.json" ]]; then
            echo "1"
        fi
    done | wc -l
}

# Main execution
main() {
    local target_group="$1"
    
    echo -e "${BLUE}AlphaFold3 MSA Data Archiving${NC}"
    echo "============================="
    echo "Base directory: $BASE_DIR"
    echo "Output directory: $OUTPUT_MSA_DIR"
    echo "Master index: $MASTER_INDEX"
    echo "Compression threads: $COMPRESSION_THREADS"
    echo
    
    if [[ $DRY_RUN -eq 1 ]]; then
        echo -e "${YELLOW}DRY RUN MODE - No archives will be created${NC}"
        echo
    fi
    
    local groups_to_process
    if [[ -n "$target_group" ]]; then
        # Process specific group
        echo "Processing specific job group: $target_group"
        groups_to_process=$(mktemp)
        echo "$target_group" > "$groups_to_process"
    else
        # Process all groups
        echo "Discovering all job groups..."
        groups_to_process=$(get_all_job_groups)
        local group_count=$(wc -l < "$groups_to_process")
        echo "Found $group_count job groups to process"
    fi
    
    echo
    
    # Initialize progress tracking
    local current_group=0
    local total_groups=$(wc -l < "$groups_to_process")
    
    # Process each group
    while read -r job_group; do
        if [[ -n "$job_group" ]]; then
            current_group=$((current_group + 1))
            echo -e "${BLUE}Processing job group [$current_group/$total_groups]: $job_group${NC}"
            
            # Find unique files for this group
            local file_list=$(find_unique_msa_files "$job_group")
            
            # Count shared MSAs for statistics
            local shared_count=$(count_shared_msas "$job_group")
            if [[ $shared_count -gt 0 ]]; then
                echo "  Jobs with shared MSAs: $shared_count"
                TOTAL_SKIPPED_FILES=$((TOTAL_SKIPPED_FILES + shared_count))
            fi
            
            # Create archive
            archive_job_group "$job_group" "$file_list"
            
            # Cleanup
            rm -f "$file_list"
            echo
        fi
    done < "$groups_to_process"
    
    # Cleanup
    rm -f "$groups_to_process"
    
    # Print summary
    echo -e "${GREEN}Archiving Summary${NC}"
    echo "================"
    echo "Archives created: $TOTAL_ARCHIVES_CREATED"
    echo "Unique MSA files archived: $TOTAL_UNIQUE_FILES"
    echo "Jobs using shared MSAs: $TOTAL_SKIPPED_FILES"
    if [[ $TOTAL_SPACE_SAVED -gt 0 ]]; then
        echo "Space saved by deduplication: $(numfmt --to=iec "$TOTAL_SPACE_SAVED")"
    fi
    echo
    
    # Query examples
    if [[ $TOTAL_ARCHIVES_CREATED -gt 0 && $DRY_RUN -eq 0 ]]; then
        echo -e "${BLUE}Query Commands for Master Index:${NC}"
        echo
        echo "# Find which archive contains a specific job's MSA:"
        echo "grep \",JOB_NAME,\" \"$MASTER_INDEX\""
        echo
        echo "# List all unique MSAs in an archive:"
        echo "grep \"^ARCHIVE_NAME\" \"$MASTER_INDEX\" | cut -d',' -f3 | sort -u"
        echo
        echo "# Count unique MSAs per job group:"
        echo "cut -d',' -f1 \"$MASTER_INDEX\" | cut -d'_' -f1 | sort | uniq -c"
    fi
}

# Script execution
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "${1:-}"
fi