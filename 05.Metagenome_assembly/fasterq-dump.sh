#!/bin/bash
#SBATCH --time=08:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=10G   # memory per CPU core
#SBATCH --job-name "fd"
#SBATCH --output "fd.log"

# load modules
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load SRA-Toolkit/3.0.3

# enable proxy to connect to internet
module load WebProxy

# Don't touch
for SAMPLE in $(cat SRA.list); do
    prefetch --max-size 100g --output-directory ./ $SAMPLE && \
	fasterq-dump --outdir ./ -t /scratch/user/kitchens/ --split-files -e 8 $SAMPLE
done
