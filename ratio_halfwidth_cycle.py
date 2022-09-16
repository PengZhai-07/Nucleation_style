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
            "#SBATCH --mail-user=zhai5108@gmail.com\n",\
            "#SBATCH --mail-type=BEGIN,END\n",\
            "#SBATCH --nodes=1\n",\
            "#SBATCH --mem=120000m"
            "#SBATCH --time=14-00:00:00\n",\
            "#SBATCH --partition=standard\n\n",\
            "#SBATCH --output=/home/%u/log/%x-%j.log\n",\
            "#SBATCH --error=/home/%u/log/error-%x-%j.log\n\n"]

par = pd.read_table("./key_par.txt", header = 0, sep=" ")
print(par.columns)
ratio = par['ratio']
halfwidth = par['halfwidth']
L = par['L']
co_size = par['co_size']
N = 8    # number of tasks for each bash script
n = 0

with open('ratio_halfwidth_cycle_part_1.sh','w') as f:       
    f.writelines(list)
    f.write("#SBATCH --job-name=gl_case1-%.0f \n" %(N))
    f.write("#SBATCH --account=yiheh0\n")
    f.write("#SBATCH --ntasks-per-node=32\n")
    f.write("# The application(s) to execute along with its input arguments and options:\n")
    f.write("# half-width(m) rigidity_ratio Lc(m)\n\n")
    for i in range(0,2):
            n = n+1
            f.write("julia --threads 4 run.jl %.6f %.6f %.6f > /home/pengzhai/log/gl_case%.0f \n\n"  \
                %(ratio.iloc[i],halfwidth.iloc[i],L.iloc[i],n))
