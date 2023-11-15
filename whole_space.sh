#!/bin/bash

#SBATCH --job-name=quartz_test
#SBATCH --array=12

##SBATCH --mail-user=pengzhai@umich.edu
##SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

#SBATCH --mem=60000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1

#SBATCH --output=log/%x-%a.log
#SBATCH --error=log/error-%x-%a.log

# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID whole_space_32.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID high_res.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID high_res_2.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID domain_size_test.txt
julia --threads ${SLURM_CPUS_PER_TASK} run.jl ${SLURM_ARRAY_TASK_ID} quartz_test.txt

# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID SSE_Creep.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID SSE_Creep_2.txt
# julia --threads 16 run.jl $SLURM_ARRAY_TASK_ID all_cases.txt