#!/bin/bash

#SBATCH -A m2_jgu-sim3		 # Account
#SBATCH -p parallel		 # Partition: parallel, smp, bigmem
#SBATCH -C skylake 		 # architecture Skylake (64 Cores) or Broadwell (40 Cores)	
#SBATCH -n 1                     # number of tasks
#SBATCH -N 1
#SBATCH --mem 200G
#SBATCH -t 10:00:00              # Run time (hh:mm:ss)


## Default Output 
WD="/lustre/miifs01/project/m2_jgu-sim3/M3-simulations/"

## Move job to Ramdisk for sufficient space
JOBDIR="/localscratch/${SLURM_JOB_ID}/"

module purge # ensures vanilla environment
module load lang/R # will load most current version of R

cp $WD/sim3.R $JOBDIR
cp -R $WD/Functions $JOBDIR
cp -R $WD/Models $JOBDIR

## Change Dir to Jobfolder
cd $JOBDIR

# Run Script
srun Rscript sim3.R -N $1 -K $2 -F $3 -R $4 -I ${SLURM_JOB_ID} -D ${WD}

