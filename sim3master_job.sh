#!/bin/bash

declare -ar nFreetime=(2)
declare -ar nRetrievals=(1000)
declare -ar OtherItems=(5)
declare -ar NPL=(16)

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

		    jobname="sim3_${N}.${K}.${nFT}.${nRet}"
      		    slurmout="Jobs/${jobname}.%j.out"
                    #echo $slurmout
		    if [[ ! -e ${slurmout} ]]; then
        		 sbatch -A "$account" -J "$jobname" -o "$slurmout" sim3_slave.sh "$N" "$K" "$nRet" "$nFT" 
        		 njobs=$((njobs + 1))

          	     fi
            done
	done
    done
done


