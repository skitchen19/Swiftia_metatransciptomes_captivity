#!/bin/bash
#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=128G   # memory per CPU core
#SBATCH --job-name "pf"
#SBATCH --output "cPF.log"

module load GCC/11.2.0
module load OpenMPI/4.1.1
module load phyloFlash/3.4.2-Python-2.7.18

DIR=/scratch/group/kitchen-group/lab/

# run phyloflash against database
for SAMPLE in $(cat sample.list.all); do
	phyloFlash.pl -lib ${SAMPLE}_rep -CPUs 16 -almosteverything -readlength 150 -taxlevel 6 \
	-read1 ${DIR}/demi/swiftia_metatranscriptomics/03_QC_filtered/${SAMPLE}_sk/${SAMPLE}_filtered_1.fastq \
	-read2 ${DIR}/demi/swiftia_metatranscriptomics/03_QC_filtered/${SAMPLE}_sk/${SAMPLE}_filtered_2.fastq
done


# run phyloflash against database

phyloFlash_compare.pl --zip MDBC01_rep.phyloFlash.tar.gz,MDBC03_rep.phyloFlash.tar.gz,MDBC04_rep.phyloFlash.tar.gz,MDBC05_rep.phyloFlash.tar.gz,MDBC06_rep.phyloFlash.tar.gz,MDBC09_rep.phyloFlash.tar.gz,MDBC11_rep.phyloFlash.tar.gz,MDBC18_rep.phyloFlash.tar.gz,MDBC21_rep.phyloFlash.tar.gz,MDBC28_rep.phyloFlash.tar.gz,MDBC29_rep.phyloFlash.tar.gz,MDBC30_rep.phyloFlash.tar.gz,MDBC32_rep.phyloFlash.tar.gz,MDBC33_rep.phyloFlash.tar.gz,MDBC34_rep.phyloFlash.tar.gz,MDBC35_rep.phyloFlash.tar.gz,MDBC36_rep.phyloFlash.tar.gz,MDBC37_rep.phyloFlash.tar.gz,MDBC38_rep.phyloFlash.tar.gz,MDBC42_rep.phyloFlash.tar.gz,MDBC44_rep.phyloFlash.tar.gz,MDBC45_rep.phyloFlash.tar.gz,MDBC46_rep.phyloFlash.tar.gz --task ntu_table --level 7 --out mergedALL
