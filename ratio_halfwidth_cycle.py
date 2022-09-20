#!/usr/bin/env python3
# -*- encoding: utf-8 -*-


# output several bash script

from re import A
import sys
import pandas as pd 

list = ["#!/bin/bash\n# The interpreter used to execute the script\n\n",\
            "# This project is to explore the phase diagram (Ru=11) with\n", \
            "# different rigidity ratio, halfwidth and characteristic slip distance \n\n",\
            "##SBATCH directives that convey submission options:\n\n",\
            "##SBATCH --mail-user=zhai5108@gmail.com\n",\
            "##SBATCH --mail-type=BEGIN,END\n",\
            "#SBATCH --nodes=1\n",\
            "#SBATCH --mem=20000m\n",\
            "#SBATCH --time=14-00:00:00\n",\
            "#SBATCH --partition=standard\n\n"]

# for high resolution, the mem should be large enough!!  default value is 20G
# but for highest resolution 20: the mem should be about 60G


par = pd.read_table("./key_par.txt", header = 0, sep=" ")
print(par.columns)
print(len(par))
ratio = par['ratio']
halfwidth = par['halfwidth']
L = par['L']
co_size = par['co_size']

n = 0              # number of tasks: calculate from the (n+1)th model
N = 8              # multi threads for each task: 4 is good!!

with open('ratio_halfwidth_cycle.sh','w') as f:       
    f.writelines(list)
    f.write("#SBATCH --account=yiheh0\n")
    f.write("#SBATCH --ntasks-per-node=%.0f \n\n" %(N))           # request N cores
    f.write("# The application(s) to execute along with its input arguments and options:\n")
    f.write("# half-width(m) rigidity_ratio Lc(m)\n\n")
    f.write("#Uncomment the following words to submit different simulations!!! \n\n\n")
    for i in range(n,len(par)):              # here from task 2
            n = n+1
            f.write("##SBATCH --job-name=glcase%.0f \n"  %(n))
            f.write("##SBATCH --output=/home/%u/log/%x-%j.log \n")
            f.write("##SBATCH --error=/home/%u/log/error-%x-%j.log \n")
            f.write("#julia --threads %.0f run.jl %.2f %.0f %.6f \n\n"\
                %(N,ratio.iloc[i],halfwidth.iloc[i],L.iloc[i]))

