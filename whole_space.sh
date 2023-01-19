#!/bin/bash

#SBATCH --array=133-210
#SBATCH --nodes=1
#SBATCH --mem=10000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case133_210
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log

julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.03
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.03
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.03
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.03
julia --threads 4 run.jl 0.8 500 0.03 4 0.00 0.03
julia --threads 4 run.jl 0.8 500 0.04 4 0.00 0.03
julia --threads 4 run.jl 0.8 500 0.05 4 0.00 0.03
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.027273
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.027273
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.027273
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.027273
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.027273
julia --threads 4 run.jl 0.8 500 0.03 4 0.00 0.027273
julia --threads 4 run.jl 0.8 500 0.04 4 0.00 0.027273
julia --threads 4 run.jl 0.8 500 0.05 4 0.00 0.027273
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.025
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.025
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.025
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.025
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.025
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.025
julia --threads 4 run.jl 0.8 500 0.03 4 0.00 0.025
julia --threads 4 run.jl 0.8 500 0.04 4 0.00 0.025
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.023077
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.023077
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.023077
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.023077
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.023077
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.023077
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.023077
julia --threads 4 run.jl 0.8 500 0.03 4 0.00 0.023077
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.021429
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.021429
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.021429
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.021429
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.021429
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.021429
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.021429
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.021429
julia --threads 4 run.jl 0.8 500 0.004 4 0.00 0.02
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.02
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.02
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.02
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.02
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.02
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.02
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.02
julia --threads 4 run.jl 0.8 500 0.003 4 0.00 0.01875
julia --threads 4 run.jl 0.8 500 0.004 4 0.00 0.01875
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.01875
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.01875
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.01875
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.01875
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.01875
julia --threads 4 run.jl 0.8 500 0.002 4 0.00 0.017647
julia --threads 4 run.jl 0.8 500 0.0025 4 0.00 0.017647
julia --threads 4 run.jl 0.8 500 0.003 4 0.00 0.017647
julia --threads 4 run.jl 0.8 500 0.004 4 0.00 0.017647
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.017647
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.017647
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.017647
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.017647
julia --threads 4 run.jl 0.8 500 0.0013 4 0.00 0.016667
julia --threads 4 run.jl 0.8 500 0.0015 4 0.00 0.016667
julia --threads 4 run.jl 0.8 500 0.002 4 0.00 0.016667
julia --threads 4 run.jl 0.8 500 0.0025 4 0.00 0.016667
julia --threads 4 run.jl 0.8 500 0.003 4 0.00 0.016667
julia --threads 4 run.jl 0.8 500 0.004 4 0.00 0.016667
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.016667
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.016667
julia --threads 4 run.jl 0.8 500 0.0006 4 0.00 0.015789
julia --threads 4 run.jl 0.8 500 0.0008 4 0.00 0.015789
julia --threads 4 run.jl 0.8 500 0.001 4 0.00 0.015789
julia --threads 4 run.jl 0.8 500 0.0013 4 0.00 0.015789
julia --threads 4 run.jl 0.8 500 0.0015 4 0.00 0.015789
julia --threads 4 run.jl 0.8 500 0.002 4 0.00 0.015789
julia --threads 4 run.jl 0.8 500 0.0025 4 0.00 0.015789
julia --threads 4 run.jl 0.8 500 0.003 4 0.00 0.015789
