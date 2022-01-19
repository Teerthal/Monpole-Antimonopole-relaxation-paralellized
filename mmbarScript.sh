#!/bin/bash
 
#SBATCH -p cluster# Send this job to the serial partition
#SBATCH -n 1                        # number of cores
#SBATCH -t 0-08:00                  # wall time (D-HH:MM)
#SBATCH -A tvachasp # Account to pull cpu hours from (commented out)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
##SBATCH --mail-user=tvachasp@asu.edu # send-to address

module load gcc/4.9.2

#make

#This file (mmbarScript.sh) is simply the "header" file that will get 
#copied out N times when the script "toRun.sh" is executed. The script 
#toRun.sh also introduces a command line (in the copied files) that 
#executes #cnSO3.out together with inputs $1, $2, $3.
#
#No need of the next 3 lines unless want to run only one set of
#paramters (so commented out):
#./cnSO3.out $1 $2 $3
#echo $1 $2 $3
#./cnSO3.out


