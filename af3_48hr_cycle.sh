#!/bin/bash

#SBATCH --job-name=af3_cycle
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --partition=normal
#SBATCH --output=af3_cycle_%A_%a.out
#SBATCH --error=af3_cycle_%A_%a.err

# Script to run AlphaFold3 submission every 2 hours for 48 hours total
# Modified to support both direct execution and SLURM submission
TOTAL_CYCLES=72
CYCLE_LENGTH=1

# Set base directory
BASE_DIR="/scratch/groups/ogozani/alphafold3"
cd "$BASE_DIR"
SUBMIT_SCRIPT="submit_dist.sh"
CYCLE_NUM=${1:-1}  # Get cycle number from command line argument, default to 1

# Check if we're running as a SLURM job or directly
if [ -z "$SLURM_JOB_ID" ]; then
    RUNNING_MODE="direct"
    echo "Running in direct mode (not as SLURM job)"
else
    RUNNING_MODE="slurm"
fi

echo "========================================="
echo "AlphaFold3 Cycle $CYCLE_NUM of $TOTAL_CYCLES"
echo "Started at $(date)"
echo "Working directory: $(pwd)"
if [ "$RUNNING_MODE" = "slurm" ]; then
    echo "SLURM Job ID: $SLURM_JOB_ID"
fi
echo "========================================="

# Check if the submit script exists (only on first run)
if [ "$CYCLE_NUM" -eq 1 ] && [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo "ERROR: Submit script not found at $SUBMIT_SCRIPT"
    echo "Please check the path and try again."
    exit 1
fi

# Make sure required scripts are executable
chmod +x "$SUBMIT_SCRIPT" submit_msa_arrays.sh batch_reuse_msa.py 2>/dev/null

# Load python module for batch_reuse_msa.py
ml python/3.12.1 2>/dev/null || module load python/3.12.1 2>/dev/null || true

# FIRST CYCLE ONLY: Run initial MSA analysis and submit array jobs
# (This is kept for compatibility when run directly via sbatch)
if [ "$CYCLE_NUM" -eq 1 ] && [ "$RUNNING_MODE" = "slurm" ]; then
    echo ""
    echo "=== First cycle: Running initial MSA analysis and submission ==="
    
    # Run batch_reuse_msa.py to identify MSA needs
    echo "Running batch_reuse_msa.py to analyze MSA requirements..."
    ./batch_reuse_msa.py
    
    if [ $? -eq 0 ]; then
        # Check if there are MSA jobs to submit
        if [ -f "msa_array_jobs.csv" ]; then
            MSA_COUNT=$(tail -n +2 msa_array_jobs.csv | grep -v "^$" | wc -l)
            
            if [ $MSA_COUNT -gt 0 ]; then
                echo ""
                echo "Found $MSA_COUNT jobs requiring MSA generation"
                echo "Submitting MSA array jobs..."
                
                ./submit_msa_arrays.sh msa_array_jobs.csv
                
                if [ $? -eq 0 ]; then
                    echo "MSA array jobs submitted successfully"
                    
                    # Check how many were actually submitted
                    MSA_RUNNING=$(squeue -u $USER -n af3_msa_array -h | wc -l)
                    echo "MSA array jobs in queue: $MSA_RUNNING"
                    
                    # If many MSA jobs were submitted, skip GPU submission this cycle
                    if [ $MSA_RUNNING -gt 100 ]; then
                        echo ""
                        echo "Large number of MSA jobs submitted ($MSA_RUNNING)"
                        echo "Skipping GPU submission this cycle to allow MSAs to start"
                        echo "GPU submissions will begin in cycle 2"
                    else
                        # Run submit_dist.sh as normal
                        echo ""
                        echo "Proceeding with GPU job submission..."
                        ./${SUBMIT_SCRIPT}
                    fi
                else
                    echo "ERROR: Failed to submit MSA array jobs"
                    echo "Continuing with GPU submission anyway..."
                    ./${SUBMIT_SCRIPT}
                fi
            else
                echo "No MSA jobs to submit (all MSAs already exist)"
                echo "Proceeding with GPU job submission..."
                ./${SUBMIT_SCRIPT}
            fi
        else
            echo "No msa_array_jobs.csv file created"
            echo "Proceeding with GPU job submission..."
            ./${SUBMIT_SCRIPT}
        fi
    else
        echo "WARNING: batch_reuse_msa.py failed"
        echo "Proceeding with GPU job submission anyway..."
        ./${SUBMIT_SCRIPT}
    fi
else
    # For cycles 2-48, just run submit_dist.sh normally
    echo "Executing: $SUBMIT_SCRIPT"
    ./${SUBMIT_SCRIPT}
fi

# Check the exit status
exit_status=$?
if [ $exit_status -eq 0 ]; then
    echo "Cycle $CYCLE_NUM completed successfully at $(date)"
else
    echo "WARNING: Cycle $CYCLE_NUM exited with status $exit_status at $(date)"
fi

# If this isn't the last cycle, submit the next job to run in $CYCLE_LENGTH hours
if [ "$CYCLE_NUM" -lt "$TOTAL_CYCLES" ]; then
    NEXT_CYCLE=$((CYCLE_NUM + 1))
    
    # Calculate start time ($CYCLE_LENGTH hours from now)
    START_TIME=$(date -d "+$CYCLE_LENGTH hours" "+%Y-%m-%dT%H:%M:%S")
    
    echo "Submitting cycle $NEXT_CYCLE to start at $START_TIME"
    
    # Submit next job with delayed start from the base directory
    cd "$BASE_DIR"
    NEXT_JOBID=$(sbatch --begin="$START_TIME" --export=ALL af3_48hr_cycle.sh $NEXT_CYCLE | awk '{print $4}')
    
    if [ $? -eq 0 ]; then
        echo "Next cycle submitted successfully with job ID: $NEXT_JOBID"
    else
        echo "ERROR: Failed to submit next cycle"
        exit 1
    fi
else
    echo "========================================="
    echo "All $TOTAL_CYCLES cycles completed!"
    echo "Final cycle finished at $(date)"
    echo "========================================="
    
    # Summary of what was processed
    if [ -f "logged_folding_jobs.csv" ]; then
        COMPLETED=$(grep -c ",0$" logged_folding_jobs.csv 2>/dev/null || echo "0")
        echo "Successfully completed jobs: $COMPLETED"
    fi
    
    # Send completion email
    if command -v mail &> /dev/null; then
        echo "AlphaFold3 48-hour cycling job completed successfully at $(date) on $(hostname)" | \
            mail -s "AlphaFold3 48hr Job Completed" dhusmann@stanford.edu
        echo "Email notification sent to dhusmann@stanford.edu"
    else
        echo "Mail command not available - please check job completion manually"
    fi
fi

echo "Cycle $CYCLE_NUM finished at $(date)"
