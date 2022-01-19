#!/bin/bash

#SBATCH -n 16                   # Use all 32 cores on the SMP node
#SBATCH -t 0-04:00                  # wall time (D-HH:MM)
##SBATCH -A tpatel28            # Account to pull cpu hours from (commented out)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=tpatel28@asu.edu # send-to address

module load intel-mpi/2018x


time srun --mpi=pmi2 "/home/tpatel28/monopoleComboPureGS/cnSO3.out" 10 10 10  

module unload intel-mpi/2018x                                                     

#cd data/

#sbatch consolidateData.sh 3000 0300

#cd ..


