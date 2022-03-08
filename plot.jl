#######################################
# Basic testing to visualize results
# #####################################

include("analyze_results.jl")   
# max slip rate
VfmaxPlot(Vfmax, t, yr2sec);
# culmulative slip
cumSlipPlot(delfsec[1:4:end,:], delfyr[1:end, :], FltX);
# 