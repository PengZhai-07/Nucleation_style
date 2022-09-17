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

#SBATCH --job-name=gl_case1-8 
#SBATCH --account=yiheh0
#SBATCH --ntasks-per-node=32

# The application(s) to execute along with its input arguments and options:
# half-width(m) rigidity_ratio Lc(m)

julia --threads 4 run.jl 0.500000 500.000000 0.010081 > /home/pengzhai/log/gl_case1 

julia --threads 4 run.jl 0.500000 1000.000000 0.011114 > /home/pengzhai/log/gl_case2 

julia --threads 4 run.jl 0.500000 2000.000000 0.011343 > /home/pengzhai/log/gl_case3 

julia --threads 4 run.jl 0.500000 5000.000000 0.011350 > /home/pengzhai/log/gl_case4 

julia --threads 4 run.jl 0.333333 500.000000 0.014247 > /home/pengzhai/log/gl_case5 

julia --threads 4 run.jl 0.333333 1000.000000 0.016496 > /home/pengzhai/log/gl_case6 

julia --threads 4 run.jl 0.333333 2000.000000 0.017008 > /home/pengzhai/log/gl_case7 

julia --threads 4 run.jl 0.333333 5000.000000 0.017025 > /home/pengzhai/log/gl_case8 

echo 'Finish submitting!'