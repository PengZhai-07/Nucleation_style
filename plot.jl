#######################################
# Basic testing to visualize results
# #####################################

using DelimitedFiles

# data storage path
turbo = "/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data"
project = "wholespace/tremor"

# the dimension of input parameters
input_parameter = readdlm("$(@__DIR__)/tremor_end_number.txt", ',', header=false)
a = size(input_parameter)[1]
# output_frequency for sliprate, stress and weakening rate
global output_freq::Int = 10   
global Domain_X::Int = 40e3

# calculate the nucleation size and plot the nucleation process
N_timestep = 600      # maximum time steps to use in sliprate to calculate nucleation size
criteria = 1e-1    # seismic threshold to measure the nucleation size
measure_threshold = 1e-3    # where measure the width of nucleation zone: 1e-7m/s for 
                            # constant weakening(expanding crack) and 1e-3m/s for fixed length patch
for index = 12
    
    # domain parameters
    Domain = input_parameter[index,1]   # amplify factor of the domain size, the current domain size is 30km*24km for 0.75 domain size
    res::Int =  input_parameter[index,2]   # resolution of mesh: should be an integer
    T::Int = input_parameter[index,3]   # total simulation time   unit:year
    # fault zone parameter
    FZlength = input_parameter[index,4]    # length of fault zone: m
    FZdepth = (Domain_X*Domain+FZlength)/2   # depth of lower boundary of damage zone  unit: m    
    halfwidth =  input_parameter[index,5]   # half width of damage zone   unit:m
    alpha = input_parameter[index,6]   # initial(background) rigidity ratio: fault zone/host rock
    cos_reduction = input_parameter[index,7]    # coseismic rigidity reduction 
    # friction parameter on fault surface
    multiple_asp::Int = input_parameter[index,8]  # effective normal stress on fault: 10MPa*multiple
    multiple_matrix = 0.1
    a_over_b = input_parameter[index,9] 
    asp_a = 0.005
    matrix_a = 0.015
    asp_b =  asp_a/a_over_b            # coseismic b increase 
    Lc= input_parameter[index,10]     # characteristic slip distance      unit:m
    matrix_asp_ratio::Int= input_parameter[index,11]     # characteristic slip distance      unit:m
    asperity_number::Int = input_parameter[index,12]  

    Fault_length = Domain*Domain_X/1000

    FILE = "$(Domain)_$(res)_$(T)_$(FZlength)_$(halfwidth)_$(alpha)_$(cos_reduction)_$(multiple_asp)_$(a_over_b)_$(Lc)_$(matrix_asp_ratio)_$(asperity_number)"
    # global FILE = "$(FZdepth)_$(halfwidth)_$(res)_$(alpha)_$(cos_reduction)_$(multiple)_$(Domain)_$(coseismic_b)_$(Lc)"
    println(FILE)                                                                                                                                                                                                                                                                                                                             
    # path to save files
    global path = "$(@__DIR__)/plots/$(project)/$(FILE)/"
    mkpath(path)
    # path to data files
    global out_path = "$(turbo)/$(project)/$(FILE)/"
    # global out_path = "$(@__DIR__)/data/$(project)/$(FILE)/"

    include("analyze_results.jl")     

    # # moment_release_example(sliprate', FltX, tStart, t, N_timestep, criteria, measure_threshold)       

    # max slip rate versus timestep
    VfmaxPlot(Vfmax, T, t)

    # healing analysis: Vfmax and regidity ratio vs. time
    healing_analysis(Vfmax, alphaa, t, yr2sec)

    # Plot friction parameters
    icsPlot(a_b, Seff, tauo, FltX,Fault_length)

    # culmulative slip
    cumSlipPlot(delfsec[1:5:end,:], delfyr[1:5:end, :], FltX, hypo, d_hypo, 1.2*T,Fault_length);
    # cumSlipPlot_no_hypocenter(delfsec[1:4:end,:], delfyr[1:end, :], FltX, 1.2*T);

    # slip rate vs timesteps
    # how many years to plot
    eqCyclePlot(sliprate', FltX, T, t,Fault_length)

    # shear stress level vs timesteps
    eqCyclePlot_stress(stress', FltX, T, t,Fault_length)
                        
    # Nucleation_example(sliprate', weakeningrate', FltX, tStart, t, N_timestep, criteria, measure_threshold)    # only plot the last seismic event

    # NS_width = Nucleation(sliprate', FltX, tStart, t, N_timestep, criteria, measure_threshold)

    # open(string(path,"nucleation info.out"), "w") do io
    #     for i = 1: size(NS_width)[1]
    #         write(io, join(NS_width[i,:], " "), "\n") 
    #     end
    # end

    


    # # # plot the variation of apparent stress
    # # # apparent_friction(stress, index_start, index_end, delfsec, index_ds_start, index_ds_end, depth, t, 2, 50)
    # apparent_friction_new(stress, index_start, index_end, delfsec_et, index_ds_start, index_ds_end, NS_width, 40, N_events)
    # apparent_friction_new_prapogation(stress, index_start, index_end, delfsec_et, index_ds_start,
    #     index_ds_end, NS_width, 40, N_events)

    # stress drop of the first artificial event 
    # # stressdrop_1(taubefore[1,:], tauafter[1,:], FltX)    # the row is the number of event

    # # coseismic stress drop for each event
    # stressdrop_2(taubefore, tauafter, FltX, tStart)    # the row is the number of event
    
    # #Plot hypo(km) vs event number, average stress drop(MPa), duration(s), moment magnitude.
    # hypo_Mw_stressdrop(hypo, Mw, del_sigma, delfafter,FltX)

    # # # sliprate versus time for the last event
    # # n = 2        # how many seismic events to plot
    # # eqCyclePlot_last_1(slipraste', FltX, tStart, t, N_timestep, n)
    # # eqCyclePlot_last_2(sliprate', FltX, tStart, t, N_timestep, n)




    # # Plot velocity dependence of b
    # velocity_dependence_b(1e-5, 1e-3, 0.019, 0.025)

    # plot the b value in rate and state friction law
    # velocity_dependence(b_value, Vfmax, t, yr2sec)
end