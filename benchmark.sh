#!/bin/bash

#SBATCH --job-name=benchmark
#SBATCH --array=13-20

##SBATCH --mail-user=pengzhai@umich.edu
##SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1

#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log

julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID benchmark.txt

# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID SSE_Creep.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID SSE_Creep_2.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID all_cases.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID whole_space_32.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID high_res.txt
