#!/bin/bash
# The interpreter used to execute the script

# This project is to explore the phase diagram (Ru=11) with
# different rigidity ratio, halfwidth and characteristic slip distance 

##SBATCH directives that convey submission options:

##SBATCH --mail-user=zhai5108@gmail.com
##SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard

#SBATCH --account=yiheh1
##SBATCH --account=lsa3
#SBATCH --ntasks-per-node=8 

# The application(s) to execute along with its input arguments and options:
# half-width(m) rigidity_ratio Lc(m) multiple(normal stress) cos_reduction

#Uncomment the following words to submit different simulations!!! 

#SBATCH --job-name=case96 
#SBATCH --output=/home/%u/log/%x-%j.log 
#SBATCH --error=/home/%u/log/error-%x-%j.log 
julia --threads 8 run.jl 0.8 500 0.025 4 0.00 0.021 


## here the 0.8 and 500 doesn't influence, because the depth of damage zone is 0 defaultly