#!/bin/bash
# The interpreter used to execute the script

# This project is to explore the phase diagram (Ru=11) with
# different rigidity ratio, halfwidth and characteristic slip distance 

##SBATCH directives that convey submission options:

##SBATCH --mail-user=zhai5108@gmail.com
##SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --mem=60000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard

#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=8 

# The application(s) to execute along with its input arguments and options:
# half-width(m) rigidity_ratio Lc(m) multiple(normal stress) cos_reduction

#Uncomment the following words to submit different simulations!!! 


#SBATCH --job-name=case92 
#SBATCH --output=/home/%u/log/%x-%j.log 
#SBATCH --error=/home/%u/log/error-%x-%j.log 
julia --threads 8 run.jl 0.8 500 0.01 4 0.00 0.019    

