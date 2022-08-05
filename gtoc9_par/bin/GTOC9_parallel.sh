#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1  
#SBATCH -n 1  
#SBATCH -c 128 
#SBATCH -J 052901
srun GTOC9_parallel