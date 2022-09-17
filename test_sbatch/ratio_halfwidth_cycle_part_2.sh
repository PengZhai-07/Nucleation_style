#!/bin/bash
# The interpreter used to execute the script

# This project is to explore the phase diagram (Ru=11) with
# different rigidity ratio, halfwidth and characteristic slip distance 

##SBATCH directives that convey submission options:

#SBATCH --mail-user=zhai5108@gmail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --mem=120000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard

#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log

#SBATCH --job-name=gl_case9-16 
#SBATCH --account=yiheh0
#SBATCH --ntasks-per-node=32

# The application(s) to execute along with its input arguments and options:
# half-width(m) rigidity_ratio Lc(m)

julia --threads 4 run.jl 0.250000 500.000000 0.018327 > /home/pengzhai/log/gl_case9 

julia --threads 4 run.jl 0.250000 1000.000000 0.021856 > /home/pengzhai/log/gl_case10 

julia --threads 4 run.jl 0.250000 2000.000000 0.022673 > /home/pengzhai/log/gl_case11 

julia --threads 4 run.jl 0.250000 5000.000000 0.022700 > /home/pengzhai/log/gl_case12 

julia --threads 4 run.jl 0.200000 500.000000 0.022366 > /home/pengzhai/log/gl_case13 

julia --threads 4 run.jl 0.200000 1000.000000 0.027206 > /home/pengzhai/log/gl_case14 

julia --threads 4 run.jl 0.200000 2000.000000 0.028338 > /home/pengzhai/log/gl_case15 

julia --threads 4 run.jl 0.200000 5000.000000 0.028375 > /home/pengzhai/log/gl_case16 

echo 'Finish submitting!'