#!/bin/bash
#SBATCH --time=08:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=64G   # memory per CPU core
#SBATCH --job-name "smr"
#SBATCH --output "cSMR.log"

DIR=/scratch/group/kitchen-group/lab/

# run sortMeRNA against database
for SAMPLE in $(cat sample.list.all); do
	${DIR}/tools/sortmerna/bin/sortmerna \
	--ref ${DIR}/tools/sortmerna/rRNA_databases_v4/bacteria.default.fasta \
	--ref ${DIR}/tools/sortmerna/rRNA_databases_v4/eukaryota.default.fasta \
	--ref ${DIR}/tools/sortmerna/rRNA_databases_v4/archaea.default.fasta \
	--ref ${DIR}/tools/sortmerna/rRNA_databases_v4/RFAM.default.fasta \
	--reads ${DIR}/demi/swiftia_metatranscriptomics/03_QC_filtered/${SAMPLE}/${SAMPLE}_filtered_1.fastq \
	--reads ${DIR}/demi/swiftia_metatranscriptomics/03_QC_filtered/${SAMPLE}/${SAMPLE}_filtered_2.fastq \
	--out2 --paired_out --fastx --aligned ${SAMPLE}_rRNA-reads --other ${SAMPLE}_non-rRNA-reads --threads 8 \
	--workdir run --blast 1
done

# create summary table
for SAMPLE in $(cat sample.list.all); do
	total=$(grep "Total reads = " ${SAMPLE}_rRNA-reads.log | awk '{print $4}')
	num=$(grep -A 4 "Coverage by" ${SAMPLE}_rRNA-reads.log |  grep -vE 'Coverage by' | awk '{print $2}' | paste -s )
	echo -e "${SAMPLE}\t$num\t$total" >> sortmerna_summary.tab
done