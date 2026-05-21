#!/bin/bash
#SBATCH --job-name=sqm_reads_test
#SBATCH --output=sqm_reads_v2.out
#SBATCH --error=sqm_reads_v2.err
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --ntasks=1

# Load modules and activate SqueezeMeta environment
module purge
module load Anaconda3/2024.02-1
source activate SqueezeMeta

# Run SqueezeMeta in sequential mode with paired reads
SqueezeMeta.pl \
  -m coassembly \
  -extassembly ./Mycoplasma_switia_exserta_bin_1_strict.fa \
  -s ./full_samples.txt \
  -f /scratch/group/kitchen-group/lab/demi/swiftia_metatranscriptomics/ \
  -p mycoplasma_expression_v2 \
  --nobins -t 8