#!/bin/bash
 
#SBATCH -J msteinhi_sims_anova   # job name
#SBATCH -p parallel              # partition: parallel, smp
#SBATCH -C skylake               # CPU type: broadwell (40 Kerne), skylake (64 Kerne)
#SBATCH -A m2_jgu-sim3           # project name
##SBATCH -N 1                    # e.g. one full node - do not do this, when your script is not using parallel code!
#SBATCH -n 1                     # Total number of tasks
#SBATCH -c 4                     # Total number of cores for the single task
#SBATCH -t 0:05:00               # Run time (hh:mm:ss)

##SBATCH --mail-user=jgoettma@uni-mainz.de
##SBATCH --mail-type=ALL

#SBATCH -o stdout_master.txt

declare -ar nFreetime(2 4)
declare -ar nRetrievals(100 250 500 100)


# now submit jobs for every permutation
# we keep track:
njobs=0

account=m2_jgu-sim3

for nFT in nFreetime $; do
   for nRet in $nRetrievals; do
      jobname="sim3${nFT}.${nRet}"
      slurmout="output/${nFT}.$nRET}.%j.out"
      echo $slurmout      
      if [[ ! -e ${slurmout} ]]; then
         sbatch -A "$account" -J "$jobname" -o "$slurmout" basic_jobscript.txt "$nRetrievals" "$nFreetime" 
         njobs=$((njobs + 1))
         if [ $njobs -gt 100 ]; then
            exit
         fi
      fi
   done
done


