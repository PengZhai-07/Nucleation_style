#################################
# Run the simulations from here
#################################

# 1. Go to par.jl and change as needed
# 2. Go to src/initialConditions/defaultInitialConditions and change as needed
# 3. Change the name of the simulation in this file
# 4. Run the simulation from terminal. (julia run.jl)
# 5. Plot results from the scripts folder

using Printf, LinearAlgebra, DelimitedFiles, SparseArrays, AlgebraicMultigrid, StaticArrays, IterativeSolvers, FEMSparse
using Base.Threads
using PyPlot
#  BLAS.set_num_threads(1)

include("$(@__DIR__)/par.jl")	    #	Set Parameters

# Put the resolution for the simulation here: should be an integer
resolution = 4

# Output directory to save data
out_dir = "$(@__DIR__)/data/immature_65_500m_no_healing/"    # healing time : 10yr
mkpath(out_dir)

P = setParameters(24e3, 500, resolution, 300)      # args = fault zone depth(m), fault zone halfwidth(m), resolution, total simulation time (years)
# println(P[3].M)

# include("$(@__DIR__)/NucleationSize.jl") 
# # calculate the nucleation size of initial rigidity ratio!!
# h_hom = NucleationSize(P)
# h_dam = h_hom/3           # with alphaa = 0.60
# CZone = CohesiveZoneSize(P)
# println("The nucleation size of homogeneous medium:", h_hom, " m")
# println("The approximate nucleation size of damage zone medium:", h_dam, " m")
# println("The Cohesive zone size of homogeneous medium:", CZone, " m")

# include("$(@__DIR__)/src/dtevol.jl")
# include("$(@__DIR__)/src/NRsearch.jl")
# include("$(@__DIR__)/src/otherFunctions.jl")

# include("$(@__DIR__)/src/main.jl")

# #simulation_time = @elapsed @time main(P, 0.65)    # all parameters, rigidity ratio 

# println("\n")

# @info("Simulation Complete!");
