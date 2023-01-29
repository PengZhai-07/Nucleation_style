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
res::Int = 100    # resolution of mesh
Domain = 0.75    # amplify factor of the domain size, the current domain size is 30km*24km for 0.75 domain size

# 4: 301 GLL nodes, average 100m on fault  
# 8: 601 GLL nodes, average 50m on fault
# 16: 1201 GLL nodes, average 25m on fault

T::Int = 50    # total simulation years 
FZdepth::Int = 0e3   # depth of lower boundary of damage zone  unit: m     20km is the maximum depth

# read the model parameters from whole_space.csv
index::Int = parse(Float64,ARGS[1])   
println(index)
# note the sequence of all imput parameters
input_parameter = readdlm("$(@__DIR__)/whole_space_3.txt", ',',  header=false)

alpha = input_parameter[index,1]   # initial(background) rigidity ratio: fault zone/host rock
halfwidth::Int =  input_parameter[index,2]   # half width of damage zone   unit:m
Lc= input_parameter[index,3]  # characteristic slip distance      unit:m
multiple::Int = input_parameter[index,4]  # effective normal stress on fault: 10MPa*multiple
cos_reduction = input_parameter[index,5]    # coseismic rigidity reduction 
coseismic_b = input_parameter[index,6]   # coseismic b increase 

println("doamin size: ",Domain)   # default is 40km*32km
println("rigidity ratio of damage zone: ",alpha)
println("halfwidth of fault zone(m): ",halfwidth)
println("characteristic slip distance(m): ", Lc)
println("effective normal stress(10MPa*multiple): ", multiple)
println("cos_reduction: ", cos_reduction)
println("cos_b: ", coseismic_b)

# vs: 2%   3%       4%         
# 0.9604   0.9409   0.9216
# 0.0396   0.0591   0.0784

turbo = "nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data"
project = "wholespace/phase_diagram_L_b"

# Output directory to save data
out_dir = "$(turbo)/$(project)/$(FZdepth)_$(halfwidth)_$(res)_$(alpha)_$(cos_reduction)_$(multiple)_$(Domain)_$(coseismic_b)_$(Lc)/"    

# clean old files 
if isdir(out_dir)
    rm(out_dir, recursive = true)
end

# To submit tens of jobs at a time(using jobarray of slurm), I need to generate the corresponding directories in advance
# See those scripts under Matlab directory
mkpath(out_dir)      

P = setParameters(FZdepth, halfwidth, res, T, alpha, multiple, Lc, Domain)   
# println(size(P[4].FltNI))   # total number of off-fault GLL nodes

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

simulation_time = @elapsed @time main(P, alpha, cos_reduction, coseismic_b)    # all parameters, rigidity ratio 

println("\n")

@info("Simulation Complete!");
