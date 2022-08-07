#!/bin/bash

declare -ar nFreetime=(4)
declare -ar nRetrievals=(1000)
declare -ar OtherItems=(4)
declare -ar NPL=(8)
declare -ar fixed_f=(0)

# submit job for every permutation
# we keep track:
njobs=0

# Guess the account
account=$(sacctmgr -n -s list user $USER format=account%30| grep -v none | head -n1 | tr -d " ")

##Spawn Jobs for each condition combination

    for N in "${OtherItems[@]}"; do
        for K in "${NPL[@]}"; do
           for nFT in "${nFreetime[@]}"; do
              for nRet in "${nRetrievals[@]}"; do
		 for fixedf in "${fixed_f[@]}"; do	

		    jobname="sim3_N${N}_K${K}_FT${nFT}_nRet${nRet}_fixedf${fixedf}"
      		    slurmout="Jobs/${jobname}.%j.out"
                    #echo $slurmout
		    if [[ ! -e ${slurmout} ]]; then
        		 sbatch -A "$account" -J "$jobname" -o "$slurmout" sim3_slave.sh "$N" "$K" "$nFT" "$nRet" "$fixedf" 
        		 njobs=$((njobs + 1))

          	     fi
	       done
            done
	done
    done
done

