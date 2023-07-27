#######################################
# Basic testing to visualize results
# #####################################

# wozhi
#FILE = "24000.0_1500_8_0.85_0.0396_3"
#FILE = "24000.0_1500_8_0.85_0.0591_5"
#FILE = "24000.0_1500_8_0.85_0.0784_7"
# Great lakes
# FILE = "20000.0_1500.0_8_0.85_0.05_4.0"   # normal stress testing
#FILE = "20000.0_500.0_20_0.85_0.05_5.0"   # resolution testing

FILE = "1.0_10_500_12000_150_0.008_5"
#FILE = "20000_500_8_0.8_0.0_4_1.0_0.023"     # Ru number testing: no healing
#FILE = "20000_500_12_0.5_0.0_5_1.0_smooth_0.23"
include("analyze_results.jl")   

# # total years to plots
N = 500

# # calculate the nucleation size and plot the nucleation process
# N_timestep = 500      # time steps to use in sliprate
# criteria = 1e-3
# NS_width = Nucleation(sliprate', FltX, tStart, t, N_timestep, criteria)

# open(string(path,"nucleation info.out"), "w") do io
#     for i = 1: size(NS_width)[1]
#         write(io, join(NS_width[i,:], " "), "\n") 
#     end
# end

# # plot the variation of apparent stress
# # apparent_friction(stress, index_start, index_end, delfsec, index_ds_start, index_ds_end, depth, t, 2, 50)
# apparent_friction_new(stress, index_start, index_end, delfsec_et, index_ds_start, index_ds_end, NS_width, 50)
# apparent_friction_new_prapogation(stress, index_start, index_end, delfsec_et, index_ds_start,
    # index_ds_end, NS_width, 40)

# max slip rate versus timestep
# VfmaxPlot(Vfmax, N, t)

# culmulative slip
cumSlipPlot(delfsec[1:4:end,:], delfyr[1:end, :], FltX, hypo, d_hypo, N);
# cumSlipPlot_no_hypocenter(delfsec[1:4:end,:], delfyr[1:end, :], FltX);

# healing analysis: Vfmax and regidity ratio vs. time
# healing_analysis(Vfmax, alphaa, t, yr2sec)

# slip rate vs timesteps
# how many years to plot
# eqCyclePlot(sliprate', FltX, N, t)





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
