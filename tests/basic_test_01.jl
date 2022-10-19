#######################################
# Basic testing to visualize results
# #####################################

using Plots
using PyPlot

# Default plot params
function plot_params()
    plt.rc("xtick", labelsize=16)
    plt.rc("ytick", labelsize=16)
    plt.rc("xtick", direction="in")
    plt.rc("ytick", direction="in")
    plt.rc("font", size=15)
    plt.rc("figure", autolayout="True")
    plt.rc("axes", titlesize=16)
    plt.rc("axes", labelsize=17)
    plt.rc("xtick.major", width=1.5)
    plt.rc("xtick.major", size=5)
    plt.rc("ytick.major", width=1.5)
    plt.rc("ytick.major", size=5)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=1.5)
    plt.rc("legend", fontsize=13)
    plt.rc("mathtext", fontset="stix")
    plt.rc("font", family="STIXGeneral")
  
    # Default width for Nature is 7.2 inches, 
    # height can be anything
    #plt.rc("figure", figsize=(7.2, 4.5))
  end

  plot_params()
  fig = PyPlot.figure(figsize=(4, 7))
  ax = fig.add_subplot(111)
  ax.hist([6,7,8,11.4], bins=collect(-0.5:16.5), orientation="horizontal", rwidth=0.9)
 
  ax.set_ylabel("Depth(km)")
  ax.set_xlabel("Event Number")
  ax.set_ylim([0,16])
  ax.invert_yaxis()
  show()
