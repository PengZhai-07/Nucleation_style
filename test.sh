#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case101
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.03

julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.03

