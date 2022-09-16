#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os 
import sys
import pandas as pd 
import numpy as np

halfwidth = [500, 1000, 2000, 5000]
alpha = [0.2, 0.25, 0.3333, 0.5]

Lc = pd.read_table("./Lc.txt", header = None)

for i in range(0,len(halfwidth)):
    for j in range(0,len(alpha)):             
        os.system("julia run.jl %.6f %.6f %.6f > log/case00.txt 2>&1 &" %(halfwidth[i],alpha[j],1e-3*Lc.iloc[i,j]))
    
# nohup julia run.jl --threads 4 > log/case01.txt 2>&1 &