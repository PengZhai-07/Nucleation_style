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


##SBATCH --job-name=case21 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 4 0.05 

##SBATCH --job-name=case22 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 4 0.06 

##SBATCH --job-name=case23 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 4 0.07 

##SBATCH --job-name=case24 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 4 0.08 

##SBATCH --job-name=case25 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 5 0.05 

##SBATCH --job-name=case26 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 5 0.06 

##SBATCH --job-name=case27 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 5 0.07 

##SBATCH --job-name=case28 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 5 0.08 

##SBATCH --job-name=case29 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 6 0.05 

##SBATCH --job-name=case30 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 6 0.06 

##SBATCH --job-name=case31 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 6 0.07 

##SBATCH --job-name=case32 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 6 0.08 

##SBATCH --job-name=case33 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 7 0.05 

##SBATCH --job-name=case34 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 7 0.06 

##SBATCH --job-name=case35 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 7 0.07 

##SBATCH --job-name=case36 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.85 1500 0.008000 7 0.08 

