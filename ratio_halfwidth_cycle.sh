#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=ratio_halfwidth_cycle
#SBATCH --mail-user=zhai5108@gmail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=14-00:00:00     
#SBATCH --account=yiheh0
#SBATCH --partition=standard


# The application(s) to execute along with its input arguments and options:

julia run.jl 500 0.2 0.008 --output=case01.out --error=case01.err    # hal-width(m) rigidity ratio Lc(m)




