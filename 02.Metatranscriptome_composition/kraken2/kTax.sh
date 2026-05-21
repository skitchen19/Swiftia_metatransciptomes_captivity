#!/bin/bash
#SBATCH --job-name=ktest
#SBATCH --output=k.out
#SBATCH --error=k.err
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --ntasks=1

# Load modules
module load GCC/13.2.0 OpenMPI/4.1.6 Kraken2/2.1.4
module load Biopython/1.84

# Change me!
DIR=/scratch/group/kitchen-group/lab/demi/swiftia_metatranscriptomics/

# run kraken2
for SAMPLE in $(cat sample.list.all); do
	kraken2 --use-names --threads 8 --db k2_db \
	--report ${SAMPLE}.report --gzip-compressed \
	--classified-out ${SAMPLE}_classified#.fastq \
    --paired ../03_QC_filtered/${SAMPLE}/${SAMPLE}_filtered_1.fastq.gz ../03_QC_filtered/${SAMPLE}/${SAMPLE}_filtered_2.fastq.gz > ${SAMPLE}.kraken
done