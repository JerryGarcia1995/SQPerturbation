#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 2
#SBATCH -n 182
#SBATCH -c 1

srun SQPerturbation.x
rm _*.dat
rm _*.txt