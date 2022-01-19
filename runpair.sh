#!/bin/bash

for twist in 00  
do
	for zm in 160
	do
    sbatch sbatchscript.sh $twist $zm 
    done
done	
