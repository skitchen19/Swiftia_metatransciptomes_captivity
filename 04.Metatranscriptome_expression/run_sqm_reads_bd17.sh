#!/bin/bash
#SBATCH --job-name=sqm_reads_test
#SBATCH --output=sqm_reads_test.out
#SBATCH --error=sqm_reads_test.err
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
  -extassembly ../bd17.fa \
  -s full_samples.txt \
  -f /scratch/group/kitchen-group/lab/demi/swiftia_metatranscriptomics/ \
  -p bd17_expression \
  --euk --nobins -t 8 --restart
