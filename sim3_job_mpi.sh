#!/bin/bash

#SBATCH -J sim3_small            # job name
#SBATCH -p parallel              # partition: parallel, smp
#SBATCH -C skylake		 # CPU type: broadwell (40 Kerne), skylake (64 Kerne)
#SBATCH -c 64
#SBATCH -A m2_jgu-sim3           # project name
#SBATCH -N 2                     # e.g. one full node - do not do this, when your script is not using parallel code!
#SBATCH -n 2                     # Total number of tasks
#SBATCH -t 5:00:00               # Run time (hh:mm:ss)

JOBDIR="/localscratch/${SLURM_JOB_ID}"

#SBATCH --mail-user=jgoettma@uni-mainz.de
#SBATCH --mail-type=ALL


module purge # ensures vanilla environment
module load lang/R # will load most current version of R

# do not forget to export OMP_NUM_THREADS, if the library you use, supports this
# not scale up to 64 threadsq
#export OMP_NUM_THREADS=64
                          
mpirun -np 2 R  --no-save -f sim3_mpi.R

