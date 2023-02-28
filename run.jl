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

global Domain_X = 40e3
global Domain_Y = 32e3

include("$(@__DIR__)/par.jl")	    #	Set Parameters

# read the model parameters from whole_space.txt
index::Int = parse(Float64,ARGS[1])   
para_file = ARGS[2]
println(index)
# note the sequence of all imput parameters
input_parameter = readdlm("$(@__DIR__)/$(para_file)", ',',  header=false)

# domain parameters
Domain::Float64 = input_parameter[index,1]   # amplify factor of the domain size, the current domain size is 30km*24km for 0.75 domain size
res::Int =  input_parameter[index,2]   # resolution of mesh: should be an integer
T::Int = input_parameter[index,3]   # total simulation time   unit:year
println("Doamin size factor: ",Domain)   # default is 40km*32km
println("Resolution: ",res)
println("Total simulation time(year): ",T)

# fault zone parameter
FZlength::Float64 = input_parameter[index,4]    # length of fault zone: m
FZdepth::Float64 = (Domain_X*Domain+FZlength)/2   # depth of lower boundary of damage zone  unit: m    
halfwidth::Float64 =  input_parameter[index,5]   # half width of damage zone   unit:m
alpha::Float64 = input_parameter[index,6]   # initial(background) rigidity ratio: fault zone/host rock
cos_reduction::Float64 = input_parameter[index,7]    # coseismic rigidity reduction 
println("Fault zone length(m): ",FZlength)   # default is 40km*32km
println("Fault zone halfwidth(m): ",halfwidth)
println("Rigidity ratio of fault zone: ",alpha)
println("Coseismic reduction of rigidity ratio: ", cos_reduction)

# friction parameter on fault surface
multiple_asp::Float64 = input_parameter[index,8]  # effective normal stress on fault: 10MPa*multiple
multiple_matrix::Float64 = 0.1              # 0.1*10MPa
a_over_b::Float64 = input_parameter[index,9] 
asp_a::Float64 = 0.005
asp_b::Float64 =  asp_a/a_over_b      # coseismic b increase 
asp_criticalness::Float64 = input_parameter[index,10]
matrix_asp_ratio::Int = input_parameter[index,11]
asperity_number::Int = input_parameter[index,12]
matrix_a::Float64 = input_parameter[index,13]

N::Int = asperity_number*(matrix_asp_ratio+1) + matrix_asp_ratio       # number of cells in RSF fault
G::Float64 = 3e10   # shear modulus of model material   unit: Pa
cell_size::Float64 = (Domain_X*Domain/2)/N
Lc::Float64 = cell_size/asp_criticalness      # nucleation size using Rice and Ruina's equation     unit:m
Dc::Float64 = pi/2*Lc/G/asp_b*(asp_b-asp_a)^2*multiple_asp*10e6     # inferred Dc value


println("Total number of cells: ", N)
println("Cell size(m): ", cell_size)
println("Effective normal stress(10MPa*multiple) in asperity: ", multiple_asp)
println("b in asperity: ", asp_b)
println("The nucleation size of homogeneous host medium(m):", Lc)
println("characteristic slip distance(m): ", Dc)
println("Cohesive zone size(m): ", 9*pi/32*G*Dc/asp_b/(multiple_asp*10e6))

# output path
turbo = "/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data"
project = "wholespace/tremor"
# Output directory to save data
out_dir = "$(turbo)/$(project)/$(Domain)_$(res)_$(T)_$(FZlength)_$(halfwidth)_$(alpha)_$(cos_reduction)_$(multiple_asp)_$(a_over_b)_$(asp_criticalness)_$(matrix_asp_ratio)_$(asperity_number)_$(matrix_a)/"    
print("Output directory: ", out_dir)
# clean old files 
if isdir(out_dir)
    rm(out_dir, recursive = true)
end
mkpath(out_dir)

P = setParameters(FZdepth::Float64, halfwidth::Float64, res::Int, T::Int, alpha::Float64, multiple_matrix::Float64,multiple_asp::Float64, Dc::Float64, Domain::Float64, asp_a::Float64, asp_b::Float64, matrix_a::Float64,matrix_asp_ratio::Int, G::Float64,N::Int,asperity_number::Int)    # usually includes constant parameters for each simulation 
# println(size(P[4].FltNI))   # total number of off-fault GLL nodes

# include("$(@__DIR__)/NucleationSize.jl") 
# # calculate the nucleation size of initial rigidity ratio!!
# h_hom_host, h_hom_dam = NucleationSize(P, alpha)
# println("The nucleation size of homogeneous host medium:", h_hom_host, " m")
# println("The nucleation size of homogeneous damage medium:", h_hom_dam, " m")
# # # h_dam = h_hom/3           # with alphaa = 0.60
# # # println("The approximate nucleation size of damage zone medium:", h_dam, " m")
# CZone = CohesiveZoneSize(P, alpha)
# println("The downlimit (damage) Cohesive zone size:", CZone, " m")

include("$(@__DIR__)/src/dtevol.jl")
include("$(@__DIR__)/src/NRsearch.jl")
include("$(@__DIR__)/src/otherFunctions.jl")

include("$(@__DIR__)/src/main.jl")

# output_frequency for sliprate, stress and weakening rate
global output_freq::Int = 10   
simulation_time = @elapsed @time main(P, alpha, cos_reduction, asp_b)     # usually includes variable parameters for each simulation 

println("\n")

@info("Simulation Complete!");
