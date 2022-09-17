#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

from re import A
import sys
import pandas as pd 
import numpy as np

list = ["##“#SBATCH” directives that convey submission options:\n\n",\
            "#SBATCH --mail-user=zhai5108@gmail.com\n",\
            "#SBATCH --mail-type=BEGIN,END\n",\
            "#SBATCH --nodes=1\n",\
            "#SBATCH --cpus-per-task=4\n",\
            "#SBATCH --time=14-00:00:00\n",\
            "#SBATCH --partition=standard\n\n",\
            "# The application(s) to execute along with its input arguments and options:\n",\
            "#SBATCH --output=/home/%u/%x-%j.log\n",\
            "#SBATCH --error=/home/%u/error-%x-%j.log\n\n",\
            "# half-width(m) rigidity_ratio Lc(m)\n\n"]

par = pd.read_table("./key_par.txt", header = 0, sep=" ")
N = 0
print(par.columns)
ratio = par['ratio']
halfwidth = par['halfwidth']
L = par['L']
co_size = par['co_size']


with open('ratio_halfwidth_cycle.sh','w') as f:     # with 关键字。优点是当子句体结束后文件会正确关闭
    f.write("#!/bin/bash\n# The interpreter used to execute the script\n\n")
    f.write("# This project is to explore the phase diagram (Ru=11) with \n\
    # different rigidity ratio, halfwidth and characteristic slip distance \n\n")
    for i in range(0,2):    
            N = N+1
            f.write("# %.0f \n" %(N))
            f.writelines(list)
            f.write("#SBATCH --job-name= gl-case%.0f \n" %(N))
            f.write("#SBATCH --account=yiheh0\n")
            f.write("julia run.jl %.6f %.6f %.6f \n\n" \
                %(ratio.iloc[i],halfwidth.iloc[i],L.iloc[i]))


