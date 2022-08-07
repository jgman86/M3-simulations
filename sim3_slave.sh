#!/bin/bash

#SBATCH -A m2_jgu-sim3		 # Account
#SBATCH -p parallel		 # Partition: parallel, smp, bigmem
#SBATCH -C skylake 		 # architecture Skylake (64 Cores) or Broadwell (40 Cores)	
#SBATCH -n 1                     # number of tasks
#SBATCH -N 1			 # allocate one full node	
#SBATCH --ramdisk=100G 		 # Reserve sufficient space for job on ramdisk	
#SBATCH -t 02:30:00              # Run time (hh:mm:ss)


## Default Output 
WD="/lustre/miifs01/project/m2_jgu-sim3/M3-simulations/"

## Move job to Ramdisk for sufficient space
JOBDIR="/localscratch/${SLURM_JOB_ID}/"
RAMDISK=$JOBDIR/ramdisk

module purge # ensures vanilla environment
module load lang/R # will load most current version of R

cp $WD/sim3.R $RAMDISK
cp -R $WD/Functions $RAMDISK
cp -R $WD/Models $RAMDISK

## Change Dir to Jobfolder
cd $RAMDISK

# Run Script
srun Rscript sim3.R -N $1 -K $2 -F $3 -R $4 -P $5 -I ${SLURM_JOB_ID} -D ${WD}

