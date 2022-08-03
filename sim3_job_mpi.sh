#!/bin/bash

#SBATCH -J sim3_mpi              # job name
#SBATCH -p parallel              # partition: parallel, smp
#SBATCH -C skylake		 # CPU type: broadwell (40 Kerne), skylake (64 Kerne)
#SBATCH -c 64
#SBATCH -A m2_jgu-sim3           # project name
#SBATCH -N 2
#SBATCH -n 2
##SBATCH --tasks-per-node=2
#SBATCH -t 18:00:00               # Run time (hh:mm:ss)

JOBDIR="/localscratch/${SLURM_JOB_ID}"

#SBATCH --mail-user=jgoettma@uni-mainz.de
#SBATCH --mail-type=ALL


module purge # ensures vanilla environment
module load lang/R # will load most current version of R
module load mpi/OpenMPI/4.1.1-GCC-11.2.0

# do not forget to export OMP_NUM_THREADS, if the library you use, supports this
# not scale up to 64 threadsq
#export OMP_NUM_THREADS=64

srun  R --no-save --slave -f  sim3_mpi.R

