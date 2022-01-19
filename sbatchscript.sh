#!/bin/bash

#SBATCH -q aggressive                   # Use all 32 cores on the SMP node
#SBATCH -N 8
#SBATCH -n 216                   # Use all 32 cores on the SMP node
#SBATCH --mem-per-cpu=4000    # Use all the RAM on the SMP node
#SBATCH -t 0-020:00                  # wall time (D-HH:MM)
##SBATCH -A drzuckerman            # Account to pull cpu hours from (commented out)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=asaurabh@asu.edu # send-to address

module load intel-mpi/2018x

#mkdir $1_$2

#cd $1_$2 

# mkdir data

time srun --mpi=pmi2 "/home/asaurabh/scratch/monopoleComboPureGS/cnSO3.out" 00 $1 $2 

module unload intel-mpi/2018x                                                     

#cd data/

#sbatch consolidateData.sh 0600

#cd ..


