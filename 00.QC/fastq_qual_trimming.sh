#!/bin/bash
#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=10G   # memory per CPU core
#SBATCH -J "fqc" # job name
#SBATCH --output "fqc"

# load modules
module load FastQC/0.11.9-Java-11

# set directories or other variables
DIR=/scratch/group/kitchen-group/lab/demi/swiftia_metatranscriptomics/

mkdir fastqc_report

#run fastqc
for SAMPLE in $(cat sample.list); do
	fastqc -o ${DIR}/03_QC_filtered/fastqc_report -f fastq -t 4 \
	${DIR}/01_raw_data/${SAMPLE}/${SAMPLE}_1.fq.gz \
	${DIR}/01_raw_data/${SAMPLE}/${SAMPLE}_2.fq.gz
done

#load modules
module purge
module load GCCcore/13.2.0
module load cutadapt/5.0

# run cutadapt v5.0
for SAMPLE in $(cat sample.list); do
        cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 50 -q 15 -j 0 \
		-o ${DIR}/03_QC_filtered/${SAMPLE}/${SAMPLE}_filtered_1.fastq \
        -p ${DIR}/03_QC_filtered/${SAMPLE}/${SAMPLE}_filtered_2.fastq \
        ${DIR}/01_raw_data/${SAMPLE}/${SAMPLE}_1.fq.gz \
        ${DIR}/01_raw_data/${SAMPLE}/${SAMPLE}_2.fq.gz
done