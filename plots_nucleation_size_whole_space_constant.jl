using DelimitedFiles
using LinearAlgebra
using StatsPlots

using Plots
using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
mpl = pyimport("matplotlib")

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
  plt.rc("font",size=12)

  # Default width for Nature is 7.2 inches, 
  # height can be anything
  #plt.rc("figure", figsize=(7.2, 4.5))

end

fig = PyPlot.figure(figsize=(10, 6));
ax = fig.add_subplot(111)

# calculate the theoretical values of nucleation size and plot the line
include("equations.jl")

gamma = pi/4
r = 1
mu = 3.2e10
mu_D = r*mu
H = 0
L=0.01
a = 0.015
sigma = 40e6
bb = range(0.017,0.032, step=0.0001)
NS_RA = zeros(length(bb))
L_b = zeros(length(bb))
for i = 1:length(bb)
    b_1 = bb[i]
    L_b[i] = 2*1.3774*mu*L/b_1/sigma
    f(h_lay) = h_lay*tanh(2*gamma*H/h_lay+ atanh(mu_D/mu)) - 2/pi*mu_D*L*b_1/sigma/(b_1-a)^2
    NS_RA[i] = Newton(f, 1000, 1e-6)
end
ax.plot(bb, NS_RA./1000, color="k")

# multiple = [4.0,5.0,6.0,7.0]
# cos_reduction = [0.05,0.06,0.07,0.08]
b = [0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.025, 0.027, 0.029, 0.031]
color_nucleation = ["r","orangered","chocolate","darkorange","orange","gold","yellow", "lawngreen", "green",
              "cyan", "dodgerblue", "blueviolet", "purple", "magenta", "pink", "black", "darkgrey", "lightgrey"]
# color can also represent the depth of hypocenters

# color_nucleation = ["r","yellow","g","b"]


for i = 1: length(b)
    N_shallow = 0
    N_deep = 0
    FILE = "0_500_8_0.8_0.0_4_0.75_$(b[i])"   # normal stress testing
    println(FILE)
    
    out_path = "$(@__DIR__)/plots/wholespace/constant b/$(FILE)/"

    NS_width = readdlm(string(out_path, "nucleation info.out"), header=false)
    # NS_width = readdlm(string(out_path, "nucleation info.out"), header=false)
    NS_average = zeros(length(b))
    NS_average[i] = mean(NS_width[2:end,2])

    # plot the theoretical nucleation size

    b_1 = b[i]
    f(h_lay) = h_lay*tanh(2*gamma*H/h_lay+ atanh(mu_D/mu)) - 2/pi*mu_D*L*b_1/sigma/(b_1-a)^2
    NS_RA_single = Newton(f, 3000, 1e-6)

    # 
    data_path = "$(@__DIR__)/data/wholespace/constant b/$(FILE)/"
    event_time = readdlm(string(data_path, "event_time.out"), header=false)
    hypo = event_time[2:end,3]
    println("Depth of all seismic events:",hypo)
    # println(size(hypo))
    # println(size(NS_width[:,1]))
    
    # # color depends on depth: shallow and deep earthquakes
    # for j = 1: length(NS_width[:,1])
    #     if hypo[j] < 30000
    #       N_shallow = N_shallow + 1
    #       if i==1&&N_shallow==1
    #         ax.plot(b[i], NS_width[j,2], "o", color = "blue",  markersize=10, label="Shallow events")
    #       else
    #         ax.plot(b[i], NS_width[j,2], "o", color = "blue",  markersize=10)
    #       end
    #     else
    #       N_deep = N_deep + 1
    #       if i==1&&N_deep==1
    #         ax.plot(b[i], NS_width[j,2], "o", color = "green", markersize=10, label="Deep events")
    #       else
    #         ax.plot(b[i], NS_width[j,2], "o", color = "green", markersize=10)
    #       end
    #     end
    # end

    # plot each earthquake for different colors
    for j = 2: length(NS_width[:,1])   # plot from the second normal seismic event(third earthquake)
      if i==1
        ax.plot(b[i], NS_width[j,2], "o", color = color_nucleation[j],  markersize=10, label = j)
      else
        ax.plot(b[i], NS_width[j,2], "o", color = color_nucleation[j],  markersize=10)    
      end    
    end

        # Plot the threotical and average values of nucleation size 
    if i==1
        ax.plot(b[i], NS_average[i], "^", color =color_nucleation[1], markerfacecolor="none", markersize=15, label = "Average value")
        ax.plot(b[i], NS_RA_single/1000, "*", color = color_nucleation[1], markerfacecolor="none",markersize=15, label = "Theoretical-R&A(2005) for a/b>0.5")
    else
        ax.plot(b[i], NS_average[i], "^", color = color_nucleation[1], markerfacecolor="none", markersize=15)
        ax.plot(b[i], NS_RA_single/1000, "*", color = color_nucleation[1], markerfacecolor="none", markersize=15)
      
    end    

end
ax.plot(bb, L_b./1000, color="k",linestyle="--", label="2*1.3774*Lb for a/b<0.3781")

ax.set_xlabel("b")
ax.set_ylabel("Nucleation size (km)")
ax.set_ylim([0,12])
ax.set_xlim([0.017,0.032])
ax.legend(loc="upper right") 
# # # healing analysis: Vfmax and regidity ratio vs. time
#healing_analysis(Vfmax, alphaa, t, yr2sec)
path = "$(@__DIR__)/plots/wholespace/constant b/"        
figname = string(path, "Nucleation size statistics.png")
fig.savefig(figname, dpi = 600, pad_inches=0.0)
show()

# event_stress = readdlm(string(out_path, "event_stress.out"), header=false)
# indx = Int(length(event_stress[1,:])/2)
# # println(length(event_stress[1,:]))  # 962
# # println(indx)          # 481
# taubefore = event_stress[:,1:indx]
# tauafter = event_stress[:,indx+1:end]

# # coseismic slip on fault for all different events(row)
# delfafter = readdlm(string(out_path, "coseismic_slip.out"), header=false)
# # print(size(delfafter))
# sliprate = readdlm(string(out_path, "sliprate.out"), header=false)
# #
# # println(size(delfafter,1))
# # println(size(delfafter,2))

# stress = readdlm(string(out_path, "stress.out"), header=false)

# start_index = get_index(stress', taubefore')
# stressdrops = taubefore .- tauafter

#         # Index of fault from 0 to 18 km
#         flt18k = findall(FltX .<= 18)[1]

#         time_vel = readdlm(string(out_path, "time_velocity.out"), header=false)
#         t = time_vel[:,1]
#         Vfmax = time_vel[:,2]
#         Vsurface = time_vel[:,3]
#         alphaa = time_vel[:,4]