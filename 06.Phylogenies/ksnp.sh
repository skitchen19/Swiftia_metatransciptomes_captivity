#!/bin/bash
#SBATCH --job-name=ks
#SBATCH --output=ks.out
#SBATCH --error=ks.err
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --ntasks=1

export PATH="/scratch/group/kitchen-group/lab/demi/swiftia_metatranscriptomics/10_metagenome/species_tree/kSNP4.1_Linux-package/kSNP4.1pkg/:$PATH"

./kSNP4.1_Linux-package/kSNP4.1pkg/MakeKSNP4infile -indir bd17_datasets -outfile ksnp_spongi.txt

./kSNP4.1_Linux-package/kSNP4.1pkg/Kchooser4 -in ksnp_spongi.txt

./kSNP4.1_Linux-package/kSNP4.1pkg/kSNP4 -k 17 -CPU 16 -outdir spongi_redo -in ksnp_spongi.txt -ML

module load GCCcore/13.2.0
module load FastTree/2.1.11

FastTree -nt -pseudo -gamma -gtr ./spongi/SNPs_all_matrix.fasta > MLtree.out
