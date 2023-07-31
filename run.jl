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

# calculate the theoretical values of nucleation size and plot the line
include("equations.jl")

include("$(@__DIR__)/par.jl")	    #	Set Parameters

# read the model parameters from whole_space.txt
index::Int = parse(Float64,ARGS[1])   
para_file = ARGS[2]
println(index)
# note the sequence of all imput parameters
input_parameter = readdlm("$(@__DIR__)/$(para_file)", ',',  header=false)

# domain parameters
Domain = input_parameter[index,1]   # amplify factor of the domain size, the current domain size is 48km*30km for 0.75 domain size
# Put the resolution for the simulation here: should be an integer
res::Int = input_parameter[index,2]     # resolution of mesh           5
# for dimain factor = 1 (48 km) 
# 4: 481 GLL nodes, average 100m on fault  
# 8: 961 GLL nodes, average 50m on fault

T::Int = input_parameter[index,3]      # total simulation years 

FZdepth::Int = input_parameter[index,4]      # depth of lower boundary of damage zone  unit: m     20km is the maximum depth
halfwidth::Int = input_parameter[index,5]   # half width of damage zone   unit:m
Lc= input_parameter[index,6]  # characteristic slip distance      unit:m
multiple::Int = input_parameter[index,7]    # effective normal stress on fault: 10MPa*multiple

cos_reduction = input_parameter[index,8]    # coseismic rigidity reduction 
coseismic_b = input_parameter[index,9]   # coseismic b increase 

# default: host rock!!
global rho1::Float64 = 2800
global vs1::Float64 = 3380.6

# # The entire medium has low rigidity
# rho1::Float64 = 2500
# vs1::Float64 = 0.6*3462
# the initial property of fualt damage zone: fault zone evolution!!!

global rho2::Float64 = 2430.6
global vs2::Float64 = 2028.36       # for healing test: define the variation of regidity in main.jl
alpha = (rho2*vs2^2)/(rho1*vs1^2)

# vs2::Float64 = sqrt(alpha)*vs1   # define the rigidity now(a constant during whole simulation)

# note: it is not necessary to define the damage zone here with healing
# But if with healing, we need to define the rigidity ratio here!!

println("domain size: ",Domain)   # default is 40km*32km
# println("rigidity ratio of damage zone: ",alpha)
println("halfwidth of fault zone(m): ",halfwidth)
println("characteristic slip distance(m): ", Lc)
println("effective normal stress(10MPa*multiple): ", multiple)
println("rigidity ratio: ", alpha)

# println("cos_reduction: ", cos_reduction)
# println("cos_b: ", coseismic_b)

# vs: 2%   3%       4%         
# 0.9604   0.9409   0.9216
# 0.0396   0.0591   0.0784

# project = "velocity_dependence_b"
project = "benchmark_project"

# Output directory to save data
out_dir = "$(@__DIR__)/data/$(project)/$(Domain)_$(res)_$(T)_$(FZdepth)_$(halfwidth)_$(Lc)_$(multiple)/"    

# # clean old files 
# if isdir(out_dir)
#     rm(out_dir, recursive = true)
# end

mkpath(out_dir)

P = setParameters(FZdepth, halfwidth, res, T, multiple, Lc, Domain)   
# println(size(P[4].FltNI))   # total number of off-fault GLL nodes

include("$(@__DIR__)/NucleationSize.jl") 
# calculate the nucleation size of initial rigidity ratio!!
h_hom_host, h_hom_dam, h_dam = NucleationSize(P, alpha)
println("The nucleation size of homogeneous host medium:", h_hom_host, " m")
println("The nucleation size of homogeneous host medium:", h_hom_dam, " m")
println("The nucleation size of layered medium with damage zone:", h_dam, " m")
# # h_dam = h_hom/3           # with alphaa = 0.60
# # println("The approximate nucleation size of damage zone medium:", h_dam, " m")
CZ_hom_host, CZ_hom_dam, CZ_dam = CohesiveZoneSize(P, alpha)
println("The Cohesive zone size of homogeneous host medium:", CZ_hom_host, " m")
println("The Cohesive zone size of homogeneous host medium:", CZ_hom_dam, " m")
println("The Cohesive zone size of layered medium with damage zone:", CZ_dam, " m")

include("$(@__DIR__)/src/dtevol.jl")
include("$(@__DIR__)/src/NRsearch.jl")
include("$(@__DIR__)/src/otherFunctions.jl")

include("$(@__DIR__)/src/main.jl")

simulation_time = @elapsed @time main(P, alpha, cos_reduction, coseismic_b)    # all parameters, rigidity ratio 

println("\n")

@info("Simulation Complete!");

