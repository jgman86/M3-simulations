#!/bin/bash

#SBATCH -A m2_jgu-sim3
#SBATCH -p parallel
#SBATCH -C broadwell
#SBATCH -n 1
##SBATCH -N 1
#SBATCH -c 16 
#SBATCH -t 10:00:00              # Run time (hh:mm:ss)


## Move job to scratch for sufficient space
JOBDIR="/localscratch/${SLURM_JOB_ID}"


module purge # ensures vanilla environment
module load lang/R # will load most current version of R

# for testing
srun Rscript sim3.R -N $1 -K $2 -F $3 -R $4 -I ${SLURM_JOB_ID} 