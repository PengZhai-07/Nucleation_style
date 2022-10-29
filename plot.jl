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

# Ru=11 test
FILE = "20000_500_4_0.5_0.0_5_1.0"     # Ru number testing: no healing

include("analyze_results.jl")   

# # max slip rate versus timestep
# VfmaxPlot(Vfmax)

# # culmulative slip
# cumSlipPlot(delfsec[1:4:end,:], delfyr[1:end, :], FltX, hypo, d_hypo);
# #cumSlipPlot_no_hypocenter(delfsec[1:4:end,:], delfyr[1:end, :], FltX);

# # # healing analysis: Vfmax and regidity ratio vs. time
# healing_analysis(Vfmax, alphaa, t, yr2sec)

# # plot the b value in rate and state friction law
# velocity_dependence(b_value, Vfmax, t, yr2sec)

# # stress drop of the first artificial event 
# # stressdrop_1(taubefore[1,:], tauafter[1,:], FltX)    # the row is the number of event

# # coseismic stress drop
# stressdrop_2(taubefore, tauafter, FltX, 6)    # the row is the number of event
# # default is the first event, so taubefore is the initial shear stress: 0.6*normal stress

# # slip rate vs timesteps
# eqCyclePlot(sliprate', FltX)


# # # sliprate versus time for the last event
# N = 500      # time steps to plot in sliprate
# n = 6           # plot the last 4 events
# Nucleation(sliprate', FltX, tStart, t, N, n)
# eqCyclePlot_last_1(sliprate', FltX, tStart, t, N, n)
# eqCyclePlot_last_2(sliprate', FltX, tStart, t, N, n)

# #Plot hypo(km) vs event number, average stress drop(MPa), duration(s), moment magnitude.
# hypo_Mw_stressdrop(hypo, Mw, del_sigma, delfafter,FltX)

# # Plot friction parameters
# icsPlot(a_b, Seff, tauo, FltX)

# Plot velocity dependence of b
velocity_dependence_b(1e-5, 1e-3, 0.019, 0.025)
