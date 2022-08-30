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
# using PyPlot    # no matplotlib in wozhi
#  BLAS.set_num_threads(1)

include("$(@__DIR__)/par.jl")	    #	Set Parameters

# Put the resolution for the simulation here: should be an integer

FZdepth = 24e3   # depth of fault zone  unit: m
halfwidth = 500   # half width of damage zone   unit:m
res = 8   # resolution of mesh
# 4: 481 GLL nodes, average 100m on fault  
# 6: 721 GLL nodes, average 67m on fault
# 8: 961 GLL nodes, average 50m on fault
# 10: 1201 GLL nodes, average 40m on fault
# 12: 1441 GLL nodes,  average 33m on fault
# 16: 1921 GLL nodes, average 25m on fault
T = 300    # total simulation years 
alpha = 0.64    # initial(background) rigidity ratio: fault zone/host rock

# vs: 2%   3%  4%              
# 0.9604   0.9409    0.9216
# 0.0396   0.0591    0.0784

# runing time
# case01: 17 hours on wozhi

cos_reduction = 0.0396    # coseismic rigidity reduction 
multiple = 3    # effective normal stress on fault: 10MPa*multiple  
# cos_reduction = 0.0591    # coseismic rigidity reduction 
# multiple = 5    # effective normal stress on fault: 10MPa*multiple
cos_reduction = 0.0784    # coseismic rigidity reduction   
multiple = 7    # effective normal stress on fault: 10MPa*multiple

# Output directory to save data
out_dir = "$(@__DIR__)/data/immature_fully_healing/$(FZdepth)_$(halfwidth)_$(res)_$(alpha)_$(cos_reduction)_$(multiple)/"    
mkpath(out_dir)

P = setParameters(FZdepth, halfwidth, res, T , alpha, multiple)   
# # println(size(P[4].FltNI))   # total number of off-fault GLL nodes

include("$(@__DIR__)/NucleationSize.jl") 
# calculate the nucleation size of initial rigidity ratio!!
h_hom_host, h_hom_dam = NucleationSize(P, alpha)
println("The nucleation size of homogeneous host medium:", h_hom_host, " m")
println("The nucleation size of homogeneous damage medium:", h_hom_dam, " m")
# # h_dam = h_hom/3           # with alphaa = 0.60
# # println("The approximate nucleation size of damage zone medium:", h_dam, " m")
CZone = CohesiveZoneSize(P, alpha)
println("The downlimit (damage) Cohesive zone size:", CZone, " m")

include("$(@__DIR__)/src/dtevol.jl")
include("$(@__DIR__)/src/NRsearch.jl")
include("$(@__DIR__)/src/otherFunctions.jl")

include("$(@__DIR__)/src/main.jl")

simulation_time = @elapsed @time main(P, alpha, cos_reduction)    # all parameters, rigidity ratio 

println("\n")

@info("Simulation Complete!");

