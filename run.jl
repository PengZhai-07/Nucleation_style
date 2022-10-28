#################################
# Run the simulations from here
#################################

# 1. Go to par.jl and change as needed
# 2. Go to src/initialConditions/defaultInitialConditions and change as needed
# 3. Change the name of the simulation in this file
# 4. Run the simulation from terminal. (julia run.jl)
# 5. Plot results from the scripts folder

using Printf, LinearAlgebra, DelimitedFiles, SparseArrays, AlgebraicMultigrid, StaticArrays, IterativeSolvers, FEMSparse
# using Distributed
using Base.Threads
# using PyPlot    # no matplotlib in wozhi
# BLAS.set_num_threads(2)  # If the underlying BLAS is using multiple threads, higher flop rates are realized

include("$(@__DIR__)/par.jl")	    #	Set Parameters

# Put the resolution for the simulation here: should be an integer

res::Int = 4   # resolution of mesh

# 4: 481 GLL nodes, average 100m on fault  
# 6: 721 GLL nodes, average 67m on fault
# 8: 961 GLL nodes, average 50m on fault
# 10: 1201 GLL nodes, average 40m on fault
# 12: 1441 GLL nodes,  average 33m on fault
# 16: 1921 GLL nodes, average 25m on fault
# 20: 1921 GLL nodes, average 20m on fault
T::Int = 200    # total simulation years 
Domain = 1.0    # amplify factor of the domain size
FZdepth::Int = 20e3   # depth of fault zone  unit: m     20km is the maximum depth

# note the sequence of all imput parameters
alpha = parse(Float64,ARGS[1])   # initial(background) rigidity ratio: fault zone/host rock
halfwidth::Int = parse(Float64, ARGS[2])   # half width of damage zone   unit:m
Lc= parse(Float64, ARGS[3])  # characteristic slip distance      unit:m
multiple::Int = parse(Float64, ARGS[4])# effective normal stress on fault: 10MPa*multiple
cos_reduction = parse(Float64, ARGS[5])    # coseismic rigidity reduction 

println("doamin size: ",Domain)   # default is 40km*32km
println("rigidity ratio of damage zone: ",alpha)
println("halfwidth of fault zone(m): ",halfwidth)
println("characteristic slip distance(m): ", Lc)
println("effective normal stress(10MPa*multiple): ", multiple)
println("cos_reduction: ", cos_reduction)

# vs: 2%   3%       4%         
# 0.9604   0.9409   0.9216
# 0.0396   0.0591   0.0784

# Output directory to save data
out_dir = "$(@__DIR__)/data/immature_fully_healing/$(FZdepth)_$(halfwidth)_$(res)_$(alpha)_$(cos_reduction)_$(multiple)_$(Domain)/"    

# clean old files 
if isdir(out_dir)
    rm(out_dir, recursive = true)
end

mkpath(out_dir)

P = setParameters(FZdepth, halfwidth, res, T, alpha, multiple, Lc, Domain)   
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

