#######################################
# Basic testing to visualize results
# #####################################
using DelimitedFiles
using LinearAlgebra

# include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")

# data storage path
global output_freq = 10
global Domain_X = 40e3


letter = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)","(l)","(m)", "(n)", "(o)", "(p)"]

turbo = "/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data"
project = "wholespace/phase_diagram_L_b"


# other input parameter
# input_parameter = readdlm("$(@__DIR__)/SSE_Creep.txt", ',',  header=false)
# input_parameter = readdlm("$(@__DIR__)/SSE_Creep_2.txt", ',',  header=false)
# input_parameter = readdlm("$(@__DIR__)/high_res.txt", ',',  header=false)
input_parameter = readdlm("$(@__DIR__)/whole_space_32_copy.txt", ',',  header=false)
a = size(input_parameter)[1]
global ii::Int = 0

plot_params()
fig = PyPlot.figure(figsize=(10, 6))

EX = [181,191,197,200].+5
N_EX = length(EX)
ax = fig.add_subplot(1, 2, 1)
ax1 = fig.add_subplot(1, 2, 2)

for index = 205                  # normal stress

    global ii += 1
    # domain parameters
    Domain = input_parameter[index,1]   # amplify factor of the domain size, the current domain size is 30km*24km for 0.75 domain size
    res::Int =  input_parameter[index,2]   # resolution of mesh: should be an integer
    T::Int = input_parameter[index,3]   # total simulation time   unit:year
    # fault zone parameter
    FZlength::Int = input_parameter[index,4]    # length of fault zone: m
    FZdepth::Int = (Domain_X*Domain+FZlength)/2   # depth of lower boundary of damage zone  unit: m    
    halfwidth::Int =  input_parameter[index,5]   # half width of damage zone   unit:m
    alpha = input_parameter[index,6]   # initial(background) rigidity ratio: fault zone/host rock
    cos_reduction = input_parameter[index,7]    # coseismic rigidity reduction 
    # friction parameter on fault surface
    multiple::Int = input_parameter[index,8]  # effective normal stress on fault: 10MPa*multiple
    a_over_b = input_parameter[index,9] 
    local a = 0.015      # a is a constant value
    coseismic_b::Float64 =  a/a_over_b            # coseismic b increase 
    Lc= input_parameter[index,10]     # characteristic slip distance      unit:m

    # global FILE = "0_500_$(res)_0.8_$(cos_reduction)_$(multiple)_$(Domain)_$(coseismic_b)_$(Lc)"
    global FILE = "$(Domain)_$(res)_$(T)_$(FZlength)_$(halfwidth)_$(alpha)_$(cos_reduction)_$(multiple)_$(a_over_b)_$(Lc)"
    println(FILE)

    global out_path = "$(turbo)/$(project)/$(FILE)/"
    # global out_path = "$(@__DIR__)/data/$(project)/$(FILE)/"

    # path to save files
    global path = "$(@__DIR__)/plots/$(project)/$(FILE)/"
        # # clean old files 
    # if isdir(path)
    #     rm(path, recursive = true)
    # end

    mkpath(path)

    include("analyze_results.jl")     

    # calculate the nucleation size and plot the nucleation process
    N_timestep::Int = input_parameter[index,11]       # time steps to use in sliprate for nucleation process
    # for i=1:6    1000

    criteria = 1e-1    # seismic threshold to measure the nucleation size
    measure_threshold = 1e-2    # where measure the width of nucleation zone: 1e-7m/s for 
                                # constant weakening(expanding crack) and 1e-3m/s for fixed length patch
        
    # theoretical nucleation size
    TNS_1 = 2/pi*mu*Lc*(a/a_over_b)/(multiple*10e6)/(a/a_over_b-a)^2/1000     # unit: km
    TNS_2 = 2*1.3774*mu*Lc/(a/a_over_b)/(multiple*10e6)/1000    # unit: km
    if a_over_b >= 0.3781
        TNS = TNS_1
    else
        TNS = TNS_2
    end

    # Nucleation_example_evolution(sliprate', weakeningrate', FltX, tStart, t, N_timestep, criteria, measure_threshold, Domain/1000, TNS, hypo/1000)    # only plot the last seismic event
    # sliprate = sliprate'
    # weakeningrate = weakeningrate'
    N = N_timestep
    Domain = Domain/1000
    hypo_1 = zeros(length(hypo), 1)
    hypo_1 = hypo./1000

    # time steps before seismic threshold 

    if 0.5 <= 5/TNS < 4  
        n_before = 160
    elseif 4 <= 5/TNS < 8 
        n_before = 150
    elseif 8 <= 5/TNS < 10
        n_before = 75
    elseif 10 <= 5/TNS < 18
        n_before = 50
    elseif 18 <= 5/TNS <= 32
        n_before = 30
    else
        n_before = 25
    end

    if ii == 3 
        n_before = 160
    elseif ii == 4
        n_before = 160
    end
    println("n_before=", n_before)

    n = length(tStart)         # how many seimsic events
    # n = 5
    # NS_width = zeros(n-1,4)
    nn = 0             # timestep when the maximum sliprate is over 1e-1m/s
    nnn::Vector{Int64} = zeros(2)
    mm = 0
    nn_trans = 0 #  when the calculation regime transforms
    nn_e = 11
    n_e = 7

    # for i = 1: n-1 
    group = 1
    NS_width = zeros(group, n_e, 4)
    
    EEP = zeros(group, n_e-1)
    EEPP = zeros(group)
    indx_1::Vector{Int64} = zeros(n_e)
    nn_before = 0
    i = 1

    #println("Time of the last seismic event(s):",tStart[end])
    indx_last = findall(t[:].<= tStart[i+5])[end]    # here i+2 is used: from the second normal event(or 3rd normal event)
    indx_last_int::Int = floor(indx_last/output_freq)
    # indx_last_int::Int = floor(indx_last)
    #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

    indx = findall(abs.(FltX) .<= Domain*Domain_X)[1]
    value = sliprate'[indx:end, indx_last_int-n_before:indx_last_int+N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
    value_1 = exp.(weakeningrate'[indx:end,indx_last_int - n_before:indx_last_int+N]).*value./1e-6    # weakening rate
    depth = FltX[indx:end]      # 10 km - 0 km
    tt = t[indx_last-n_before*output_freq:10:indx_last+N*output_freq]./yr2sec

    indx_around_hypocenter = findall((hypo_1[i+5] .- TNS) .<= abs.(depth) .<= (hypo_1[i+5] .+ TNS))  # get the indx around nucleation location

    # println("Index around hypocenter: ", indx_around_hypocenter)
    # find the index of maximum speed is over 1e-1m/s
    for j = 1: n_before     
        if maximum(value[indx_around_hypocenter,j]) >= 1e-8
            mm = j
            break
        end
    end

    println("The number of the timestep over 1e-8 m/s is:", mm)

    for k = 1:n_e

        # find the timestep when sliprate' first exceeds 1e-1 m/s, numbered from 1e-8 m/s
        for j = 2:N
            if  maximum(value[indx_around_hypocenter,j]) >= criteria/(10^(n_e-k))
                nn = j
                break
            end
        end
        println("The number of the timestep from 1e-8 m/s when maximum sliprate is over ",criteria/(10^(n_e-k)) ,"m/s is:", nn - mm)

        # measure the width of nucleation zone

        indx_nucleation = findall(value[indx_around_hypocenter, nn] .>= criteria/(10^(n_e-k+1)))       # using the second line(n_before+2) to define the width of nucleation size
        # println(indx_nucleation)
        new_depth = new_depth = depth[indx_around_hypocenter][indx_nucleation]      # here the depth should be this!!
        downdip_depth = maximum(new_depth)
        updip_depth = minimum(new_depth)

        NS_width[i,k, 1] = (downdip_depth + updip_depth)/2       # middle point 
        NS_width[i,k, 2] = downdip_depth - updip_depth        # width of nucleation zone
        NS_width[i,k, 3] = updip_depth                   # updip
        NS_width[i,k, 4] = downdip_depth          # downdip

        indx_1[k] = findall(FltX .>= NS_width[i,k,1])[end]   

    end
    # the center of the asperity
    println("indx_1=", indx_1)

    for j = 2:N
        if  value[indx_1[end],j] >= criteria/(10^(nn_e-1))
            nnn[1] = j
            break
        end
    end

    for j = 2:N
        if  value[indx_1[end],j] >= criteria
            nnn[2] = j
            break
        end
    end
    println("The index:", nnn)

    ax.plot(tt[nnn[2]].-tt[nnn[1]:nnn[2]], value[indx_1[end], nnn[1]:nnn[2]], label= string("a/b=",round(a_over_b, digits=4), "&\$RD_{RS}=\$", round(1000*Lc/5, digits=4)) , color="k") 
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("Time to instability(year)")
    ax.set_ylabel("V(m/s)")
    ax.legend(loc="lower left")
    ax.set_title("(b)", loc="left")

    # ax1.plot(tt[nnn[2]].-tt[nnn[1]:nnn[2]], value_1[indx_1[end], nnn[1]:nnn[2]], label= string("a/b=",round(a_over_b, digits=4), "&\$RD_{RS}=\$", round(1000*Lc/5, digits=4)) ) 
    # ax1.set_yscale("log")
    # ax1.set_xscale("log")
    # ax1.set_xlabel("Time to instability(year)")
    # ax1.set_ylabel("Ω")
    # # ax1.legend(loc="upper right")
    # ax1.set_title("(b)", loc="left")

end

path = "$(@__DIR__)/plots/wholespace/phase_diagram_L_b"        
figname = string(path, "_foreshock_representative_nucleation_style_continuous_example.png")
fig.savefig(figname, dpi = 600)
#show()
