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

FZdepth = 24e3   # depth of fault zone
halfwidth = 1000
res = 6   # resolution of mesh
# 4: 481 GLL nodes, average 100m on fault  
# 6: 721 GLL nodes, average 67m on fault
# 8: 961 GLL nodes, average 50m on fault
# 10: 1201 GLL nodes, average 40m on fault
# 12: 1441 GLL nodes,  average 33m on fault
# 16: 1921 GLL nodes, average 25m on fault
T = 300    # total simulation years 
alpha = 0.85    # velocity ratio: fault zone/host rock

# Output directory to save data
out_dir = "$(@__DIR__)/data/fully_healing/$(FZdepth)_$(halfwidth)_$(res)_$(alpha)/"    
mkpath(out_dir)

P = setParameters(FZdepth, halfwidth, res, T , alpha)   
# println(size(P[4].FltNI))

include("$(@__DIR__)/NucleationSize.jl") 
# calculate the nucleation size of initial rigidity ratio!!
h_hom = NucleationSize(P)
println("The nucleation size of homogeneous medium:", h_hom, " m")
# h_dam = h_hom/3           # with alphaa = 0.60
# println("The approximate nucleation size of damage zone medium:", h_dam, " m")
CZone = CohesiveZoneSize(P, alpha)
println("The Cohesive zone size:", CZone, " m")

include("$(@__DIR__)/src/dtevol.jl")
include("$(@__DIR__)/src/NRsearch.jl")
include("$(@__DIR__)/src/otherFunctions.jl")

include("$(@__DIR__)/src/main.jl")

simulation_time = @elapsed @time main(P, alpha)    # all parameters, rigidity ratio 

println("\n")

@info("Simulation Complete!");

