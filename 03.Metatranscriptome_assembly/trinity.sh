#!/bin/bash
#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=260G   # memory per CPU core
#SBATCH -J "build"   # job name
#SBATCH --output "build2.log"

#load modules
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load Trinity/2.15.1

dir=/scratch/group/kitchen-group/lab/demi/swiftia_metatranscriptomics/06_transcriptome_assembly

# command to run
Trinity --seqType fq  \
--left $dir/ALL_non-rRNA-reads_fwd.fq \
--right $dir/ALL_non-rRNA-reads_rev.fq \
--CPU 16 --max_memory 200G --min_contig_length 300 \
--min_kmer_cov 2

#load modules
module purge
module load GCC/13.2.0
module load CD-HIT/4.8.1

# command to run
cd-hit-est -i allReads.Trinity.fasta -o Swiftia_cdhit.fa -c 0.95 -n 10 -d 0 -M 0 -T 8