#######################################
# Basic testing to visualize results
# #####################################

using DelimitedFiles

# Universal parameter
res = 16    # resolution of mesh
Domain = 0.75    # amplify factor of the domain size, the current domain size is 30km*24km for 0.75 domain size
T = 500    # total simulation years 
FZdepth = 0   # depth of lower boundary of damage zone  unit: m     20km is the maximum depth

# other input parameter
input_parameter = readdlm("$(@__DIR__)/whole_space_1.txt", ',',  header=false)
a = size(input_parameter)[1]

for index = 59:78
    alpha = input_parameter[index,1]   # initial(background) rigidity ratio: fault zone/host rock
    halfwidth::Int =  input_parameter[index,2]   # half width of damage zone   unit:m
    Lc= input_parameter[index,3]  # characteristic slip distance      unit:m
    multiple::Int = input_parameter[index,4]  # effective normal stress on fault: 10MPa*multiple
    cos_reduction = input_parameter[index,5]    # coseismic rigidity reduction 
    coseismic_b = input_parameter[index,6]   # coseismic b increase 

    global FILE = "$(FZdepth)_$(halfwidth)_$(res)_$(alpha)_$(cos_reduction)_$(multiple)_$(Domain)_$(coseismic_b)_$(Lc)"
    println(FILE)
    include("analyze_results.jl")     

    # # total years to plotss
    N = 300          

    # calculate the nucleation size and plot the nucleation process
    N_timestep = 2000      # time steps to use in sliprate
    criteria = 1e-1    # seismic threshold to measure the nucleation size
    measure_threshold = 1e-3    # where measure the width of nucleation zone: 1e-7m/s for 
                                # constant weakening(expanding crack) and 1e-3m/s for fixed length patch

    # moment_release_example(sliprate', FltX, tStart, t, N_timestep, criteria, measure_threshold)       
                        
    # Nucleation_example(sliprate', FltX, tStart, t, N_timestep, criteria, measure_threshold)    # only plot the last seismic event

    # NS_width = Nucleation(sliprate', FltX, tStart, t, N_timestep, criteria, measure_threshold)

    # open(string(path,"nucleation info.out"), "w") do io
    #     for i = 1: size(NS_width)[1]
    #         write(io, join(NS_width[i,:], " "), "\n") 
    #     end
    # end

    # max slip rate versus timestep
    VfmaxPlot(Vfmax, N, t)

    # culmulative slip
    #cumSlipPlot(delfsec[1:4:end,:], delfyr[1:end, :], FltX, hypo, d_hypo, N);
    cumSlipPlot_no_hypocenter(delfsec[1:4:end,:], delfyr[1:end, :], FltX, N);

    # healing analysis: Vfmax and regidity ratio vs. time
    healing_analysis(Vfmax, alphaa, t, yr2sec)

    # slip rate vs timesteps
    # how many years to plot
    eqCyclePlot(sliprate', FltX, N, t)


    # # # plot the variation of apparent stress
    # # # apparent_friction(stress, index_start, index_end, delfsec, index_ds_start, index_ds_end, depth, t, 2, 50)
    # apparent_friction_new(stress, index_start, index_end, delfsec_et, index_ds_start, index_ds_end, NS_width, 40, N_events)
    # apparent_friction_new_prapogation(stress, index_start, index_end, delfsec_et, index_ds_start,
    #     index_ds_end, NS_width, 40, N_events)

    # stress drop of the first artificial event 
    # # stressdrop_1(taubefore[1,:], tauafter[1,:], FltX)    # the row is the number of event

    # # coseismic stress drop
    # stressdrop_2(taubefore, tauafter, FltX, tStart)    # the row is the number of event
    # # default is the first event, so taubefore is the initial shear stress: 0.6*normal stress

    # #Plot hypo(km) vs event number, average stress drop(MPa), duration(s), moment magnitude.
    # hypo_Mw_stressdrop(hypo, Mw, del_sigma, delfafter,FltX)

    # # # sliprate versus time for the last event
    # # n = 2        # how many seismic events to plot
    # # eqCyclePlot_last_1(sliprate', FltX, tStart, t, N_timestep, n)
    # # eqCyclePlot_last_2(sliprate', FltX, tStart, t, N_timestep, n)


    # # Plot friction parameters
    # icsPlot(a_b, Seff, tauo, FltX)

    # # Plot velocity dependence of b
    # velocity_dependence_b(1e-5, 1e-3, 0.019, 0.025)

    # # plot the b value in rate and state friction law
    # velocity_dependence(b_value, Vfmax, t, yr2sec)
end