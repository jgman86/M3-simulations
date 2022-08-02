#!/bin/bash

#SBATCH -A m2_jgu-sim3
#SBATCH -p parallel
#SBATCH -C broadwell
#SBATCH -n 1                    
#SBATCH -c 1
#SBATCH -t 00:00:20              # Run time (hh:mm:ss)


module purge # ensures vanilla environment
module load lang/R # will load most current version of R

# for testing
srun Rscript --vanilla --no-save --slave  sim_test.R  -N $1 -K $2  -F $3  -R $4 
