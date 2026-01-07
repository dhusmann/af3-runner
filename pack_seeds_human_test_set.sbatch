#!/bin/bash
#SBATCH --job-name=pack-seeds-hts
#SBATCH --output=logs/pack-seeds_hts_%j.out
#SBATCH --error=logs/pack-seeds_hts_%j.err
#SBATCH --partition=hns
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail
cd /scratch/groups/ogozani/alphafold3
THREADS=2 PARALLEL_DIRS=16 ./compress_seeds.sh output/human_test_set
