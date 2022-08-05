#!/bin/bash

declare -ar nFreetime=(2 4)
declare -ar nRetrievals=(100 250 500 1000)
declare -ar OtherItems=(1 2 3 4 5)
declare -ar NPL=(1 2 4 8 16)
declare -ar scale_b=(0.1)

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
		 for baseline in "${scale_b[@]}"; do	

		    jobname="sim3_N${N}_K${K}_FT_${nFT}_nRet${nRet}_scale_$baseline"
      		    slurmout="Jobs/${jobname}.%j.out"
                    #echo $slurmout
		    if [[ ! -e ${slurmout} ]]; then
        		 sbatch -A "$account" -J "$jobname" -o "$slurmout" sim3_slave.sh "$N" "$K" "$nFT" "$nRet" "$baseline" 
        		 njobs=$((njobs + 1))

          	     fi
	       done
            done
	done
    done
done


