#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case101
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.03

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case102
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.03

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case103
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.03

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case104
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.03

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case105
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.03 4 0.00 0.03

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case106
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.04 4 0.00 0.03

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case107
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.05 4 0.00 0.03

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case108
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.027273

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case109
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.027273

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case110
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.027273

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case111
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.027273

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case112
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.027273

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case113
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.03 4 0.00 0.027273

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case114
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.04 4 0.00 0.027273

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case115
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.05 4 0.00 0.027273

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case116
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.025

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case117
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.025

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case118
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.025

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case119
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.025

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case120
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.025

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case121
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.025

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case122
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.03 4 0.00 0.025

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case123
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.04 4 0.00 0.025

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case124
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.023077

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case125
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.023077

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case126
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.023077

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case127
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.023077

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case128
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.023077

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case129
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.023077

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case130
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.023077

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case131
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.03 4 0.00 0.023077

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case132
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.021429

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case133
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.021429

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case134
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.021429

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case135
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.021429

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case136
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.021429

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case137
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.021429

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case138
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.021429

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case139
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.025 4 0.00 0.021429

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case140
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.004 4 0.00 0.02

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case141
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.02

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case142
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.02

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case143
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.02

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case144
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.02

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case145
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.02

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case146
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.016 4 0.00 0.02

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case147
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.02 4 0.00 0.02

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case148
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.003 4 0.00 0.01875

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case149
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.004 4 0.00 0.01875

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case150
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.01875

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case151
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.01875

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case152
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.01875

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case153
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.01875

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case154
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.012 4 0.00 0.01875

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case155
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.002 4 0.00 0.017647

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case156
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0025 4 0.00 0.017647

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case157
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.003 4 0.00 0.017647

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case158
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.004 4 0.00 0.017647

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case159
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.017647

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case160
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.017647

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case161
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.008 4 0.00 0.017647

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case162
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.01 4 0.00 0.017647

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case163
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0013 4 0.00 0.016667

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case164
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0015 4 0.00 0.016667

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case165
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.002 4 0.00 0.016667

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case166
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0025 4 0.00 0.016667

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case167
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.003 4 0.00 0.016667

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case168
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.004 4 0.00 0.016667

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case169
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.005 4 0.00 0.016667

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case170
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.006 4 0.00 0.016667

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case171
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0006 4 0.00 0.015789

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case172
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0008 4 0.00 0.015789

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case173
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.001 4 0.00 0.015789

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case174
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0013 4 0.00 0.015789

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case175
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0015 4 0.00 0.015789

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case176
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.002 4 0.00 0.015789

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case177
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.0025 4 0.00 0.015789

#SBATCH --nodes=1
#SBATCH --mem=20000m
#SBATCH --time=14-00:00:00
#SBATCH --partition=standard
#SBATCH --account=yiheh1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=case178
#SBATCH --output=/home/%u/log/%x-%j.log
#SBATCH --error=/home/%u/log/error-%x-%j.log
julia --threads 4 run.jl 0.8 500 0.003 4 0.00 0.015789

