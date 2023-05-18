using DelimitedFiles
using LinearAlgebra
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

using DelimitedFiles

include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")


# data storage path
turbo = "/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data"
project = "wholespace/tremor"
FILE = ["0.125_32_1_0.0_0.0_1.0_0.0_0.3_0.5_0.4_2_10_0.015","0.125_32_1.0_5000.0_500.0_0.7225_0.0_0.25_0.5_0.4_2_10_0.015_0.001"]
println(FILE)
# Global variables
yr2sec = 365*24*60*60

xx = zeros(length(FILE), 10000::Integer)
VV = zeros(length(FILE), 10000::Integer)
indx::Vector{Int64} = zeros(length(FILE))
max_V::Vector{Float64} = zeros(length(FILE))
min_V::Vector{Float64} = zeros(length(FILE))
mean_V::Vector{Float64} = zeros(length(FILE))

for i = 1:length(FILE)

  out_path = "$(turbo)/$(project)/$(FILE[i])/"
  # Read data
  time_vel = readdlm(string(out_path, "time_velocity.out"), header=false, Float64)
  t = time_vel[:,1]         # all real timsteps
  Vfmax = time_vel[:,2]
  alphaa = time_vel[:,3]         # initial background rigidity ratio
  #b_value = time_vel[:,5]

  # Order of storage: Seff, tauo, FltX, cca, ccb, xLf
  params = readdlm(string(out_path, "params.out"), header=false)
  Seff = params[1,:]
  tauo = params[2,:]
  FltX = params[3,:]
  #println("Dimension of FltX:",size(FltX))
  cca = params[4,:]
  ccb = params[5,:]
  a_b = params[6,:]
  Lc = params[7,:]

  event_time = readdlm(string(out_path, "event_time.out"), header=false)
  tStart = event_time[:,1]
  println("Start time of all seismic events(s):",tStart) 
  tEnd = event_time[:,2]
  println("Duration of all seismic events(s):",tEnd-tStart)

  hypo = event_time[:,3]
  d_hypo = event_time[:,4]    # unit: m 
  println("Cumulative slips when earthquakes happen:",d_hypo) 
  println("Depth of all seismic events:",hypo)

  x = t./yr2sec*12*30
  indx_S = 11
  indx_E = 15
  indx_1::Int64 = findall(t .<= tStart[indx_S])[end]
  indx_2::Int64 = findall(t .<= tEnd[indx_E])[end]
  # indx_first = x[indx_1:indx_2][1]
  # indx_last = x[indx_1:indx_2][end]
  indx[i]= indx_2 - indx_1 + 1
  println(indx[i])
  xx[i,1:indx[i]] = x[indx_1:indx_2].- x[indx_1]
  VV[i,1:indx[i]] = Vfmax[indx_1:indx_2]
  max_V[i] = maximum(VV[i,1:indx[i]]) 
  min_V[i] = minimum(VV[i,1:indx[i]])
  mean_V[i] = mean(VV[i,1:indx[i]])
end

plot_params()
fig = PyPlot.figure(figsize=(10, 6));
ax = fig.add_subplot(111)

time_length_1 = xx[1,indx[1]] - xx[1,1]
println(time_length_1)

for i = 1:length(FILE)

  ax.plot(xx[i, 1:indx[i]].*(time_length_1/(xx[i,indx[i]] .- xx[i,1])), (log.(VV[i,1:indx[i]]).- log(min_V[i]))/(log(max_V[i])- log(min_V[i])), lw = 2.0)
  ax.set_ylim([0, 1])
  # ax.plot(xx[i, 1:indx[i]].*(time_length_1/(xx[i,indx[i]] .- xx[i,1])), VV[i,1:indx[i]], lw = 2.0, label= i)
  # ax.set_yscale("log")

end

# ax.plot([indx_first, indx_last],[1e-1, 1e-1] , "k", linestyle="-", label="Regular earthquake threshold")
# ax.plot([indx_first, indx_last],[Vthres, Vthres], "k", linestyle="--", label="Inertial iterm threshold")
# ax.plot([indx_first, indx_last],[1e-7, 1e-7], "g", linestyle=":", label="Tremor threshold")

ax.set_xlabel("Time (day)")
ax.set_ylabel("Max. Slip rate (m/s)")
ax.set_title("Normalized STF")

# #  ax.legend([lab1, lab2], loc=0)
show()

figname = string("Normalized STF.png")
fig.savefig(figname, dpi = 300)