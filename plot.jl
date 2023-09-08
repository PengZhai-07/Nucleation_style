#######################################
# Basic testing to visualize results
# #####################################

using DelimitedFiles

# data storage path
global output_freq = 10
global Domain_X = 40e3

turbo = "/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data"
project = "wholespace/phase_diagram_L_b/"

# other input parameter
# input_parameter = readdlm("$(@__DIR__)/SSE_Creep.txt", ',',  header=false)     # seismic events: 1-12,16,25,28,31,36,39,42
# input_parameter = readdlm("$(@__DIR__)/SSE_Creep_2.txt", ',',  header=false)
# input_parameter = readdlm("$(@__DIR__)/high_res.txt", ',',  header=false)
input_parameter = readdlm("$(@__DIR__)/whole_space_32.txt", ',',  header=false)
# input_parameter = readdlm("$(@__DIR__)/high_res_2.txt", ',',  header=false)
# input_parameter = readdlm("$(@__DIR__)/domain_size_test.txt", ',',  header=false)

a = size(input_parameter)[1]

for index = 131

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

    # N = 500
    N = T

    # calculate the nucleation size and plot the nucleation process
    N_timestep::Int = input_parameter[index,11]       # time steps to use in sliprate for nucleation process
    # for i=1:6    1000

    criteria = 1e-1    # seismic threshold to measure the nucleation size
    measure_threshold = 1e-2    # where measure the width of nucleation zone: 1e-7m/s for 
                                # constant weakening(expanding crack) and 1e-2 m/s for fixed length patch
      

    # moment_release_example(sliprate', FltX, tStart, t, N_timestep, criteria, measure_threshold, Domain/1000)       

    # sliprate_analysis(Vfmax, alphaa, t, yr2sec, sliprate', FltX, N, Domain/1000)


     # Plot friction parameters
    icsPlot(a_b, Seff, tauo, FltX.-Domain*Domain_X/1000/2, Domain/1000)

    # # max slip rate versus timestep
    # VfmaxPlot(Vfmax, N, t)

    # # healing analysis: Vfmax and regidity ratio vs. time
    # V_max = healing_analysis(Vfmax, alphaa, t, yr2sec)
    # open(string(path,"Maximum_slip_velocity.out"), "w") do io
    #     write(io, join(V_max), "\n") 
    # end

    # # slip rate vs timesteps
    # # how many years to plot
    # eqCyclePlot(sliprate', FltX, N, t, Domain/1000)
   
    
    # # culmulative slip
    # cumSlipPlot(delfsec[1:4:end,:], delfyr[1:4:end, :], FltX, hypo, d_hypo, N, Domain/1000);
    # # cumSlipPlot_no_hypocenter(delfsec[1:4:end,:], delfyr[1:end, :], FltX, N);

    # # theoretical nucleation size
    # if a_over_b >= 0.3781
    #     TNS = 2/pi*mu*Lc*(a/a_over_b)/(multiple*10e6)/(a/a_over_b-a)^2/1000     # unit: km
    # else
    #     TNS = 2*1.3774*mu*Lc/(a/a_over_b)/(multiple*10e6)/1000    # unit: km
    # end
    # println("Theoretical nucleation size: ", TNS)

    # # Nucleation_example(sliprate', weakeningrate', FltX, tStart, t, N_timestep, criteria, measure_threshold, Domain/1000)    # only plot the last seismic event
    
    # EEP = Nucleation_example_evolution(sliprate', weakeningrate', FltX, tStart, t, N_timestep, criteria, measure_threshold, Domain/1000, TNS, hypo/1000, [3,5,6,8,9])    # only plot the last seismic event
    # # EEP = Nucleation_example_evolution(sliprate', weakeningrate', FltX, tStart, t, N_timestep, criteria, measure_threshold, Domain/1000, TNS, hypo/1000, characteristic_index[3:7])    # only plot the last seismic event

    # ## Nucleation_example_no_weakening_rate(sliprate', FltX, tStart, t, N_timestep, criteria, measure_threshold)    # only plot the last seismic event

    # NS_width, min_Ω = Nucleation(sliprate', weakeningrate', FltX, tStart, t, N_timestep, criteria, measure_threshold, Domain/1000)
    # println(NS_width)
    # println(min_Ω)

    # df = mean(abs.(NS_width[:,1].-(Domain_X/1e3*Domain/2)))/(Domain_X/1e3*Domain/4)     # deviation factor of hypocenter
    # println(df)
    # if  0 <= df < 0.2
    #     rupture_style = "Symmetric-bilateral rupture"
    # elseif 0.2 <= df <= 0.6
    #     rupture_style = "Asymmetric-bilateral rupture"
    # elseif 0.6 < df 
    #     rupture_style = "Unilateral rupture"
    # end
    # println(rupture_style)
    # Ω = mean(min_Ω)
    # if  EEP < 0         # this criterion is not accurate
    #     nucleation_style = "fixed length nucleation"
    # else
    #     nucleation_style = "constant weakening nucleation"
    # end
    # println(nucleation_style)

    # open(string(path,"nucleation.out"), "w") do io
    #     for i = 1: size(NS_width)[1]
    #         write(io, join(vcat(NS_width[i,:],min_Ω[i])," "), "\n") 
    #     end
    # end
    # open(string(path,"nucleation_stats.out"), "w") do io
    #     write(io, join(df), "\n") 
    #     write(io, join(rupture_style), "\n") 
    #     write(io, join(Ω), "\n") 
    #     write(io, join(EEP), "\n")
    #     write(io, join(nucleation_style), "\n") 
    # end

    # # sliprate versus time for the last event
    # n = 2        # how many seismic events to plot
    # eqCyclePlot_last_1(sliprate', FltX, tStart, t, floor(Int64, N_timestep*1.2), n, Domain/1000, characteristic_index[2:n+1], TNS, vs2/1000)
    # eqCyclePlot_last_2(sliprate', FltX, tStart, t, floor(Int64, N_timestep*1.2), n, Domain/1000, characteristic_index[2:n+1], TNS, vs2/1000)
    # eqCyclePlot_last_3(sliprate', FltX, tStart, t, floor(Int64, N_timestep*1.2), n, Domain/1000, characteristic_index[2:n+1], TNS, vs2/1000)
    

    #  coseismic stress drop
    # stressdrop_2(taubefore, tauafter, FltX, tStart, Domain/1000)    # the row is the number of event
    # # default is the first event, so taubefore is the initial shear stress: 0.6*normal stress
    
    # # # plot the variation of apparent stress
    # # # apparent_friction(stress, index_start, index_end, delfsec, index_ds_start, index_ds_end, depth, t, 2, 50)
    # apparent_friction_new(stress, index_start, index_end, delfsec_et, index_ds_start, index_ds_end, NS_width, 40, N_events)
    # apparent_friction_new_prapogation(stress, index_start, index_end, delfsec_et, index_ds_start,
    #     index_ds_end, NS_width, 40, N_events)

    # stress drop of the first artificial event 

    # # stressdrop_1(taubefore[1,:], tauafter[1,:], FltX)    # the row is the number of event

    # #Plot hypo(km) vs event number, average stress drop(MPa), duration(s), moment magnitude.
    # hypo_Mw_stressdrop(hypo, Mw, del_sigma, delfafter,FltX)



    # # Plot velocity dependence of b
    # velocity_dependence_b(1e-5, 1e-3, 0.019, 0.025)

    # # plot the b value in rate and state friction law
    # velocity_dependence(b_value, Vfmax, t, yr2sec)
    
    mpl.pyplot.close()

end

