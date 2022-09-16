#!/bin/bash
# The interpreter used to execute the script

# This project is to explore the phase diagram (Ru=11) with 
    # different rigidity ratio, halfwidth and characteristic slip distance 

# 1 
##“#SBATCH” directives that convey submission options:

#SBATCH --mail-user=zhai5108@gmail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard

# The application(s) to execute along with its input arguments and options:
#SBATCH --output=/home/%u/%x-%j.log
#SBATCH --error=/home/%u/error-%x-%j.log

# half-width(m) rigidity_ratio Lc(m)

#SBATCH --job-name= gl_case1 
#SBATCH --account=yiheh0
julia run.jl 0.500000 500.000000 0.010081 

# 2 
##“#SBATCH” directives that convey submission options:

#SBATCH --mail-user=zhai5108@gmail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard

# The application(s) to execute along with its input arguments and options:
#SBATCH --output=/home/%u/%x-%j.log
#SBATCH --error=/home/%u/error-%x-%j.log

# half-width(m) rigidity_ratio Lc(m)

#SBATCH --job-name= gl_case2 
#SBATCH --account=yiheh0
julia run.jl 0.500000 1000.000000 0.011114 

