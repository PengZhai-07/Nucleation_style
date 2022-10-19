#!/bin/bash
# The interpreter used to execute the script

# This project is to explore the phase diagram (Ru=11) with
# different rigidity ratio, halfwidth and characteristic slip distance 

##SBATCH directives that convey submission options:

##SBATCH --mail-user=zhai5108@gmail.com
##SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --mem=30000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard

#SBATCH --account=yiheh0
#SBATCH --ntasks-per-node=8 

# The application(s) to execute along with its input arguments and options:
# half-width(m) rigidity_ratio Lc(m) multiple(normal stress) cos_reduction

#Uncomment the following words to submit different simulations!!! 


##SBATCH --job-name=case41
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.5 100 0.0070 5 0.00 

##SBATCH --job-name=case42 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.5 200 0.0081 5 0.00 

#SBATCH --job-name=case47 
#SBATCH --output=/home/%u/log/%x-%j.log 
#SBATCH --error=/home/%u/log/error-%x-%j.log 
julia --threads 8 run.jl 0.5 500 0.010081 5 0.00  

##SBATCH --job-name=case44 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.5 1000 0.011114 5 0.00 

##SBATCH --job-name=case45 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.5 2000 0.011343  5 0.00 

##SBATCH --job-name=case46 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.5 5000 0.011350  5 0.00 
