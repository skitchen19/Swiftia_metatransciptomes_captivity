#!/bin/bash
#SBATCH --job-name=sqm
#SBATCH --output=sqm.out
#SBATCH --error=sqm.err
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
  -s genome_samples.txt \
  -f ./ \
  -p genomeSM \
  --euk -t 8 --cleaning
