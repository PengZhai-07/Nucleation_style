#################################
# Run the simulations from here
#################################

# 1. Go to par.jl and change as needed
# 2. Go to src/initialConditions/defaultInitialConditions and change as needed
# 3. Change the name of the simulation in this file
# 4. Run the simulation from terminal. (julia run.jl)
# 5. Plot results from the scripts folder

using Printf, LinearAlgebra, DelimitedFiles, SparseArrays, AlgebraicMultigrid, StaticArrays, IterativeSolvers, FEMSparse
using BenchmarkTools
import AlgebraicMultigrid as AMG

# using MUMPS     # for MPI on multi cores
# using Distributed
# using Octavian
using Base.Threads
# using PyPlot    # no matplotlib in wozhi
BLAS.set_num_threads(8)  # If the underlying BLAS(!) is using multiple threads, higher flop rates are realized
# functions in LinearAlgebra with BLAS end with !

global Domain_X = 40e3
global Domain_Y = 32e3
global output_freq = 10
global tol::Float64 = 1e-6      
# 1e-6 is good enough, even though for a smaller tol(1e-14), CG can keep converging!

println("tol: ", tol)

# output path
turbo = "/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data"
# need to change this for different project 
project = "wholespace/velocity_ground"


include("$(@__DIR__)/par.jl")	    #	Set Parameters

# read the model parameters from whole_space.txt
index::Int = parse(Float64,ARGS[1])   
para_file = ARGS[2]
println(index)
# note the sequence of all imput parameters
input_parameter = readdlm("$(@__DIR__)/$(para_file)", ',',  header=false)

# domain parameters
Domain = input_parameter[index,1]   # amplify factor of the domain size, the current domain size is 30km*24km for 0.75 domain size
res::Int =  input_parameter[index,2]   # resolution of mesh: should be an integer
T::Int = input_parameter[index,3]   # total simulation time   unit:year
println("Doamin size factor: ",Domain)   # default is 40km*32km
println("Resolution: ",res)
println("Total simulation time(year): ",T)

# fault zone parameter
FZlength::Int = input_parameter[index,4]    # length of fault zone: m
FZdepth::Int = (Domain_X*Domain+FZlength)/2   # depth of lower boundary of damage zone  unit: m    
halfwidth::Int =  input_parameter[index,5]   # half width of damage zone   unit:m
alpha = input_parameter[index,6]   # initial(background) rigidity ratio: fault zone/host rock
cos_reduction = input_parameter[index,7]    # coseismic rigidity reduction 
println("Fault zone length(m): ",FZlength)   # default is 40km*32km
println("Fault zone halfwidth(m): ",halfwidth)
println("Rigidity ratio of fault zone: ",alpha)
println("Coseismic reduction of rigidity ratio: ", cos_reduction)

# friction parameter on fault surface
multiple::Int = input_parameter[index,8]  # effective normal stress on fault: 10MPa*multiple
a_over_b = input_parameter[index,9] 
coseismic_a = 0.015
coseismic_b =  coseismic_a/a_over_b            # coseismic b increase 
Lc = input_parameter[index,10]     # characteristic slip distance      unit:m
IDstate::Int = input_parameter[index,12]     # evolution of state variable: 2: aging law   3: slip law
println("Effective normal stress(10MPa*multiple): ", multiple)
println("Coseismic b: ", coseismic_b)
println("characteristic slip distance(m): ", Lc)
println("Evolution of state variable(2:aging, 3: slip):",IDstate) 

# Output directory to save data
out_dir = "$(turbo)/$(project)/$(Domain)_$(res)_$(T)_$(FZlength)_$(halfwidth)_$(alpha)_$(cos_reduction)_$(multiple)_$(a_over_b)_$(Lc)_$(IDstate)/"    
println("Output directory: ", out_dir)
# clean old files 
if isdir(out_dir)
    rm(out_dir, recursive = true)
end

mkpath(out_dir)      

P = setParameters(FZdepth, halfwidth, res, T, alpha, multiple, Lc, Domain, IDstate, coseismic_a, coseismic_b)    # usually includes constant parameters for each simulation 
# println(size(P[4].FltNI))   # total number of off-fault GLL nodes

include("$(@__DIR__)/NucleationSize.jl")
# calculate the nucleation size of initial rigidity ratio!!
h_hom_host, h_hom_dam = NucleationSize(P, alpha, coseismic_a, coseismic_b)
println("The nucleation size of homogeneous host medium:", h_hom_host, " m")
println("The nucleation size of homogeneous damage medium:", h_hom_dam, " m")
# # h_dam = h_hom/3           # with alphaa = 0.60
# # println("The approximate nucleation size of damage zone medium:", h_dam, " m")
CZone = CohesiveZoneSize(P, alpha, coseismic_b)
println("The downlimit (damage) Cohesive zone size:", CZone, " m")

include("$(@__DIR__)/src/dtevol.jl")
include("$(@__DIR__)/src/NRsearch.jl")
include("$(@__DIR__)/src/otherFunctions.jl")

include("$(@__DIR__)/src/main.jl")

simulation_time = @elapsed @time main(P, alpha, cos_reduction, coseismic_b, Domain)     # usually includes variable parameters for each simulation 
println("simulation_time(hours): ", simulation_time/3600)
println("simulation_time(days): ", simulation_time/3600/24)

println("\n")

@info("Simulation Complete!");
