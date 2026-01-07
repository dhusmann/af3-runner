#!/bin/bash
#SBATCH -J pack-seeds
#SBATCH -c 16
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -o logs/pack-seeds.%j.out
#SBATCH -e logs/pack-seeds.%j.err
#SBATCH -p hns

# Ensure env on compute node
# source ~/.bashrc
# conda activate base 2>/dev/null || true
# module load pigz 2>/dev/null || true

export THREADS=${SLURM_CPUS_PER_TASK:-12}
export COMPRESSION_LEVEL=7
export PARALLEL_DIRS=1     # pigz uses all threads for one dir at a time
export COMPRESSOR=auto

#./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/XXX
echo compressing /scratch/groups/ogozani/alphafold3/output/yeast/EFM
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/yeast/EFM
echo compressin ./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/yeast/RKM
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/yeast/RKM
echo compressing ./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/yeast/SET
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/yeast/SET
echo compressing /scratch/groups/ogozani/alphafold3/output/human_EF_KMTs/EF2
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/human_EF_KMTs/EF2
echo compressing /scratch/groups/ogozani/alphafold3/output/human_EF_KMTs/EEF1A
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/human_EF_KMTs/EEF1A
echo compressing /scratch/groups/ogozani/alphafold3/output/ecoli/ecRPL11
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/ecoli/ecRPL11
echo compressing /scratch/groups/ogozani/alphafold3/output/ecoli/ecEFTU
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/ecoli/ecEFTU
echo compressing /scratch/groups/ogozani/alphafold3/output/SETD6
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/SETD6
echo compressing /scratch/groups/ogozani/alphafold3/output/RPL36A
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/RPL36A
echo compressing /scratch/groups/ogozani/alphafold3/output/RPL29
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/RPL29
echo compressing /scratch/groups/ogozani/alphafold3/output/NSD2i
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/NSD2i
echo compressing /scratch/groups/ogozani/alphafold3/output/dpEEF1A
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/dpEEF1A
echo compressing /scratch/groups/ogozani/alphafold3/output/arabidopsis_EEF1A
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/arabidopsis_EEF1A
echo compressing /scratch/groups/ogozani/alphafold3/output/human_test_set
./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/human_test_set
