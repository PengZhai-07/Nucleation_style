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

#SBATCH --account=yiheh0
#SBATCH --ntasks-per-node=8 

# The application(s) to execute along with its input arguments and options:
# half-width(m) rigidity_ratio Lc(m)

#Uncomment the following words to submit different simulations!!! 


##SBATCH --job-name=glcase1 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.500000 500.000000 0.010081 

##SBATCH --job-name=glcase2 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.500000 1000.000000 0.011114 

##SBATCH --job-name=glcase3 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.500000 2000.000000 0.011343 

##SBATCH --job-name=glcase4 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.500000 5000.000000 0.011350 

##SBATCH --job-name=glcase5 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.333333 500.000000 0.014247 

##SBATCH --job-name=glcase6 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.333333 1000.000000 0.016496 

##SBATCH --job-name=glcase7 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.333333 2000.000000 0.017008 

##SBATCH --job-name=glcase8 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.333333 5000.000000 0.017025 

##SBATCH --job-name=glcase9 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.250000 500.000000 0.018327 

##SBATCH --job-name=glcase10 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.250000 1000.000000 0.021856 

##SBATCH --job-name=glcase11 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.250000 2000.000000 0.022673 

##SBATCH --job-name=glcase12 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.250000 5000.000000 0.022700 

##SBATCH --job-name=glcase13 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.200000 500.000000 0.022366 

##SBATCH --job-name=glcase14 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.200000 1000.000000 0.027206 

##SBATCH --job-name=glcase15 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.200000 2000.000000 0.028338 

##SBATCH --job-name=glcase16 
##SBATCH --output=/home/%u/log/%x-%j.log 
##SBATCH --error=/home/%u/log/error-%x-%j.log 
#julia --threads 8 run.jl 0.200000 5000.000000 0.028375 

