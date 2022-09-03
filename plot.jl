#######################################
# Basic testing to visualize results
# #####################################
FILE = "24000.0_500_8_0.85_0.0396_3"
# FILE = "24000.0_500_4_0.64_0.0591_5"
# FILE = "24000.0_500_8_0.64_0.0784_7"

include("analyze_results.jl")   

# max slip rate
VfmaxPlot(Vfmax, t, yr2sec);

# culmulative slip
#cumSlipPlot(delfsec[1:4:end,:], delfyr[1:end, :], FltX, hypo, d_hypo);
cumSlipPlot(delfsec[1:end,:], delfyr[1:end, :], FltX);
# cumSlipPlot(delfsec[1:4:end,:], delfyr[1:end, :], FltX, hypo);

# healing analysis: Vfmax and regidity ratio vs. time
healing_analysis(Vfmax, alphaa, t, yr2sec)

# coseismic stress drop
stressdrop_2(taubefore[1,:], tauafter[1,:], FltX)    # the row is the number of event
# default is the first event, so taubefore is the initial shear stress: 0.6*normal stress

# Plot slip vs event number
#slipPlot(delfafter', rupture_len, FltX, Mw, tStart)           # delfafter: 179*481

# slip rate vs timesteps
eqCyclePlot(sliprate', FltX)

# plot alphaa vs time
# alphaaPlot(alphaa, t, yr2sec)

# Plot friction parameters
icsPlot(a_b, Seff, tauo, FltX)