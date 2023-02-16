#!/bin/bash

#SBATCH --job-name=tremor_test
#SBATCH --array=5-6

##SBATCH --mail-user=pengzhai@umich.edu
##SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20

#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1

#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log

julia --threads 20 run.jl $SLURM_ARRAY_TASK_ID tremor_end_number.txt
