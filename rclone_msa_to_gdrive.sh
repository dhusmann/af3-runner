#!/bin/bash

#SBATCH --job-name=RclMSA
#SBATCH --output=logs/Rclone_msa_%j.out
#SBATCH --error=logs/Rclone_msa_%j.err
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=24:00:00

ml system
ml rclone/1.59.1

rclone copy /scratch/groups/ogozani/alphafold3/jobs gozani_labshare_alphafold:/alphafold3/msa/ \
  --include "*data.json" \
  --include "*msa.json" \
  --progress \
  --transfers 5 \
  --checkers 2 \
  --fast-list 
