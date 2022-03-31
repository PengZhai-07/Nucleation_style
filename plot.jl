#######################################
# Basic testing to visualize results
# #####################################

include("analyze_results.jl")   

# max slip rate
# VfmaxPlot(Vfmax, t, yr2sec);

# culmulative slip
# cumSlipPlot(delfsec[1:4:end,:], delfyr[1:end, :], FltX);

# healing analysis
# healing_analysis(Vfmax, alphaa, t, yr2sec)

# coseismic stress drop
#stressdrop_2(taubefore[84,:], tauafter[84,:], FltX)    # the row is the number of event

#slipPlot(delfafter', rupture_len, FltX, Mw, tStart)           # delfafter: 179*481

# slip rate vs timesteps
#eqCyclePlot(sliprate', FltX)

# plot alphaa
#alphaaPlot(alphaa, t, yr2sec)

# Plot friction parameters
icsPlot(a_b, Seff, tauo, FltX)