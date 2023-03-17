#!/bin/bash

#SBATCH --job-name=Phase_diagram
#SBATCH --array=1

##SBATCH --mail-user=39 pengzhai@umich.edu
##SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

#SBATCH --mem=5000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1

#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log

julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID whole_space.txt

