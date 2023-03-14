##############################
#  PLOTTING SCRIPTS
##############################

using Plots
using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
# mpl = @pyimport("matplotlib")

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

# Plot friction parameters
function icsPlot(a_b, Seff, tauo, FltX,Fault_length)
    plot_params()
    fig = PyPlot.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)            
    
    ax.plot(Seff, FltX, "b", label="Normal Stress")
    ax.plot(tauo, FltX, "orange",linestyle="-", label="Shear Stress")
    ax.set_xlabel("Stresses (MPa)")
    ax.set_ylabel("Depth (km)")
    ax.legend(loc="lower right") 
    
    col="tab:green"
    ax2 = ax.twiny()
    ax2.plot(a_b, FltX, "g",label="(a-b)")
    #println(FltX)
    # seismogenic_depth = findall(abs(10) .< abs.(FltX) .<= abs(20))   # note: unit of FltX here is km
    # a_b[seismogenic_depth] .= a_b[seismogenic_depth] .- 0.006
    # #println(a_b)
    # ax2.plot(a_b, FltX, "r",label="coseismic (a-b)", linestyle=":")
    ax2.set_xlabel("Rate-state friction value (a-b)", color=col)
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)  
    ax2.set_xlim([-0.01,0.030])
    ax2.legend(loc="lower left") 

    ax.set_ylim([0,Fault_length])
    ax.invert_yaxis()
    show()
    
    figname = string(path, "initial_condition.png")
    fig.savefig(figname, dpi = 300)
end

# Plot Vfmax
function VfmaxPlot(Vfmax, N, t)
    plot_params()
    fig = PyPlot.figure(figsize=(15, 5))
    ax = fig.add_subplot(211)

    t_seconds = (N) * 365 * 24 * 60 * 60 
    indx_last = floor(Int, findall(t .<= t_seconds)[end])   # last event!

    ax.plot(Vfmax, lw = 2.0)
    ax.plot([0, indx_last],[1e-1, 1e-1] , "k", linestyle="-", label="Regular earthquake threshold")
    ax.plot([0, indx_last],[1e-3, 1e-3] , "k", linestyle="--", label="Inertial iterm threshold")
    ax.plot([0, indx_last],[1e-7, 1e-7] , "k", linestyle=":", label="Tremor threshold")
    ax.legend(loc="upper right") 
    ax.set_xlabel("Time steps")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.set_ylim([1e-10,1e2])
    # show()
    
    ax = fig.add_subplot(212)
    Vfmax_ex = Vfmax[4000:6000]
    ax.plot(Vfmax_ex, lw = 2.0)
    ax.plot([0, length(Vfmax_ex)],[1e-1, 1e-1] , "k", linestyle="-", label="Regular earthquake threshold")
    ax.plot([0, length(Vfmax_ex)],[1e-3, 1e-3] , "k", linestyle="--", label="Inertial iterm threshold")
    ax.plot([0, length(Vfmax_ex)],[1e-7, 1e-7] , "k", linestyle=":", label="Tremor threshold")
    ax.legend(loc="upper right") 
    ax.set_xlabel("Time steps")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.set_ylim([1e-10,1e2])

    figname = string(path, "Vfmax.png")
    fig.savefig(figname, dpi = 300)
end

# Plot alpha and Vfmax on the same plot
function healing_analysis(Vf, alphaa, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    x = t./yr2sec
    indx_last = x[end]

    ax.plot(x, Vf, lw = 2.0)
    ax.plot([0, indx_last],[1e-1, 1e-1] , "k", linestyle="-", label="Regular earthquake threshold")
    ax.plot([0, indx_last],[1e-3, 1e-3] , "k", linestyle="--", label="Inertial iterm threshold")
    ax.plot([0, indx_last],[1e-6, 1e-6] , "k", linestyle=":", label="Tremor threshold")
    ax.legend(loc="upper right") 
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    #ax.set_xlim([0, 600])
    ax.set_ylim([1e-10, 1e2])
    
    col="tab:red"
    ax2 = ax.twinx()
    
    ax2.plot(t./yr2sec, alphaa.*100, c=col, lw=2.0, label="Shear modulus ratio")   
    ax2.set_ylabel("Shear Modulus (% of host rock)")
    ax2.set_ylim([20, 110])
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)

    #  ax.legend([lab1, lab2], loc=0)
    # show()
    
    figname = string(path, "healing_analysis.png")
    fig.savefig(figname, dpi = 300)
end

# Plot cumulative slip
function cumSlipPlot(delfsec, delfyr, FltX, hypo, d_hypo, N,Fault_length)
    
    indx_1 = findall(abs.(FltX) .<= Fault_length*7/8)[1]
    indx_2 = findall(abs.(FltX) .>= Fault_length/8)[end]

    print(indx_1)
    print(indx_2)

    delfsec2 = transpose(delfsec[:, indx_1:indx_2])
    delfyr2 = transpose(delfyr)

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))

    ax = fig.add_subplot(111)
    plt.rc("font",size=12)
    ax.plot(delfyr2, FltX, color="royalblue", lw=1.0)
    ax.plot(delfsec2, FltX[indx_1:indx_2], color="chocolate", lw=1.0)
    ax.plot(d_hypo, hypo./1000, "*", color="saddlebrown", markersize=10)
    ax.set_xlabel("Cumulative Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,Fault_length])

    L = (N) * 365 * 24 * 60 * 60 * 1e-9
    ax.set_xlim([0,L])
    #ax.set_xlim([0,9.0])
    ax.invert_yaxis()

    # show()
    
    figname = string(path, "cumulative_slip.png")
    fig.savefig(figname, dpi = 300)

end

# sliprate versus time plot
# N is the number of years to plot
function eqCyclePlot(sliprate, FltX, N, t,Fault_length)

    t_seconds = N * 365 * 24 * 60 * 60 
    indx_last = findall(t .<= t_seconds)[end]   # last event!
    indx_last_int::Int = floor(indx_last/output_freq)

    indx = findall(abs.(FltX) .<= 30)[1]
    value = sliprate[indx:end,1:indx_last_int]
    
    # print(findall(abs.(sliprate) .<= 1e-12))
    value[findall(abs.(value) .<= 1e-15)] .= 1e-15 

    # depth = FltX[indx:end]

    plot_params()
    fig = PyPlot.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    c = ax.imshow(value, cmap="rainbow", aspect="auto",
                  norm=matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e-6),
                  interpolation="none",    # the interpolation method decide the final slip rate distrbution!!
                  extent=[0,length(value[1,:]), 0,Fault_length])

    # for stress
    #  c = ax.imshow(value, cmap="inferno", aspect="auto",
                  #  vmin=22.5, vmax=40,
                  #  interpolation="bicubic",
                  #  extent=[0,length(seismic_slipvel[1,:]), 0,16])
    
    ax.set_xlabel("Variable Timesteps")
    ax.set_ylabel("Depth (km)")

    ax.invert_yaxis()
    cbar = fig.colorbar(c, label = "Slip rate(m/s)")
    #   cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    # show()
    figname = string(path, "sliprate_time_1.png")
    fig.savefig(figname, dpi = 600)

    plot_params()
    fig = PyPlot.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    c = ax.imshow(value, cmap="rainbow", aspect="auto",
                  norm=matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e-5),
                  interpolation="none",    # the interpolation method decide the final slip rate distrbution!!
                  extent=[0,length(value[1,:]), 0,Fault_length])

    # for stress
    #  c = ax.imshow(value, cmap="inferno", aspect="auto",
                  #  vmin=22.5, vmax=40,
                  #  interpolation="bicubic",
                  #  extent=[0,length(seismic_slipvel[1,:]), 0,16])
    
    ax.set_xlabel("Variable Timesteps")
    ax.set_ylabel("Depth (km)")

    ax.invert_yaxis()
    cbar = fig.colorbar(c, label = "Slip rate(m/s)")
    #   cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    # show()
    figname = string(path, "sliprate_time_2.png")
    fig.savefig(figname, dpi = 600)
    
end

# Plot the stressdrops after each earthquake
function stressdrop_2(taubefore, tauafter, FltX, tStart, Fault_length, multiple_asp)
    N = length(tStart)
    plot_params()
    fig = PyPlot.figure(figsize=(15, 30));
    n = 8
    for i = 2:2+n-1 
    ax = fig.add_subplot(n/2,2,i-1)
      ax.plot(taubefore[i,:], FltX, lw = 2.0, color="tab:orange", 
              label="Shear stress before the earthquake", alpha=1.0);
      ax.plot(tauafter[i,:], FltX, lw = 2.0, color="tab:blue", 
              label="Shear stress after the earthquake", alpha=1.0);
      ax.set_xlabel("Stress drop (MPa)");
      ax.set_ylabel("Depth (km)");
      ax.set_ylim([0,Fault_length]);
      ax.set_xlim([0.55*multiple_asp*10, 0.65*multiple_asp*10]);
      ax.invert_yaxis();
      plt.legend();
    end
      show()
      figname = string(path, "shear_stress_following.png")
      fig.savefig(figname, dpi = 600)
    #  end
end

function eqCyclePlot_stress(stress, FltX, N, t,Fault_length, multiple_asp)

    t_seconds = N * 365 * 24 * 60 * 60 
    indx_last = findall(t .<= t_seconds)[end]   # last event!
    indx_last_int::Int = floor(indx_last/output_freq)

    indx = findall(abs.(FltX) .<= 30)[1]
    value = stress[indx:end,1:indx_last_int]    

    # depth = FltX[indx:end]

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)

    c = ax.imshow(value, cmap="viridis", aspect="auto",
                norm=matplotlib.colors.Normalize(0.60*multiple_asp*10,maximum(maximum(value))),
                interpolation="none",    # the interpolation method decide the final slip rate distrbution!!
                extent=[0,length(value[1,:]), 0,Fault_length])

    # for stress
    #  c = ax.imshow(value, cmap="inferno", aspect="auto",
                  #  vmin=22.5, vmax=40,
                  #  interpolation="bicubic",
                  #  extent=[0,length(seismic_slipvel[1,:]), 0,16])
    
    ax.set_xlabel("Variable Timesteps")
    ax.set_ylabel("Depth (km)")

    ax.invert_yaxis()
    cbar = fig.colorbar(c, label = "Shear stress(MPa)")
    #   cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    # show()
    figname = string(path, "stress_time.png")
    fig.savefig(figname, dpi = 600)
    
end
 
function moment_release_example(sliprate, FltX, tStart, t, N, criteria, measure_threshold)
    n_before = 150             # 200 time steps before seismic threshold 
    n = length(tStart)         # how many seimsic events
    # n = 5
    nn = 0             # timestep when the maximum sliprate is over 1e-1m/s
    plot_params()
    fig = PyPlot.figure(figsize=(10, 8))
    # for i = 1: n-1 
    
    for i = n-3        # plot the i_th normal earthquake: choose by yourself    if i=n-1, then plot the last one!!
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:] .<= tStart[i+1])[end]   
        indx_last_int::Int = floor(indx_last/output_freq)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)
        t_coseismic = t[indx_last:output_freq:indx_last+N*output_freq]   # accurate time for each sliprate
        
        indx = findall(abs.(FltX) .<= 30)[1]
        depth = FltX[indx:end]
        value = sliprate[indx:end, indx_last_int:indx_last_int+N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
        
        crack_length = zeros(length(t_coseismic))
        V = zeros(length(t_coseismic))
        moment_rate = zeros(length(t_coseismic))
        moment_acceleration = zeros(length(t_coseismic))
        rupture_speed = zeros(length(t_coseismic))

        # find the timestep when sliprate first exceeds 1e-1 m/s
        for j = 2:N
            if maximum(value[:, j]) >= criteria       # 1e-1 m/s
                nn = j
                println("The maximum sliprate is over 0.1m/s")
                break
            end
        end
        t_coseismic = t_coseismic .- t_coseismic[nn]

        moment_before = 0
        # calculate the crack length
        for k in 2:length(t_coseismic)
            # measure the width of nucleation zone for each timestep
            indx_nucleation = findall(value[:, k] .>= measure_threshold)       # using the second line(n_before+2) to define the width of nucleation size
            new_depth = FltX[indx:end][indx_nucleation]
            downdip_depth = maximum(new_depth)
            updip_depth = minimum(new_depth)
            V[k] = maximum(value[:, k])
            # expanding crack
            crack_length[k] = downdip_depth - updip_depth        # width of nucleation zone
             # moment release: assuming that the rupture width is the same with rupture length
            moment_rate[k] = 3.2e10 * (crack_length[k]*1000)^2 * mean(value[indx_nucleation, k])   # unit: Nm/s
            if k < nn
                moment_before = moment_before + 0.5*(moment_rate[k] + moment_rate[k-1])*(t_coseismic[k] - t_coseismic[k-1])
            end
        end
        for k in eachindex(t_coseismic)[2:end-1]
            rupture_speed[k] = (crack_length[k+1] - crack_length[k-1])/(t_coseismic[k+1] - t_coseismic[k-1])   # central difference   
            moment_acceleration[k] = (moment_rate[k+1] - moment_rate[k-1])/(t_coseismic[k+1] - t_coseismic[k-1])   # central difference
        end
        # plot slip rate profile
        # ax = fig.add_subplot(n-1, 1, i)
        ax = fig.add_subplot(211)
        # println(size(t[indx_last_int:indx_last_int + N]))
        # println(size(value))
        
        ax.plot(t_coseismic[2:end], crack_length[2:end], color="red", label="crack length")        # plot every five steps
        ax.plot(0, crack_length[nn], "*",color="black", markersize=15, label="measured nucleation size")
        ax.set_xlabel("t(s)")
        ax.set_ylim([0, 15])
        ax.set_ylabel("Crack length(km)")
        ax.legend(loc="upper left")

        col="tab:blue"
        ax2 = ax.twinx()
        # index_middle::Int = floor(size(sliprate)[1]/2)+1
        # println(index_middle)
        # ax2.plot(t_coseismic[2:end], value[index_middle, 2:end], color="green", label="slip rate at the middle of nucleation zone")    
        # ax2.plot(t_coseismic[2:end], V[2:end], color="green", label="Maximum slip rate")    
        # ax2.get_xaxis().set_tick_params(color=col)
        # ax2.tick_params(axis="x", labelcolor=col) 
        # ax2.set_yscale("log")
        # ax2.set_ylim([1e-4,1e0])    
        # # ax2.set_ylabel("Slip rate at the middle of nucleation zone(m/s)")
        # ax2.set_ylabel("Maximum slip rate(m/s)")
        # ax2.legend(loc="upper right")

        ax2.plot(t_coseismic[2:end], rupture_speed[2:end], "o",color=col, label="rupture speed", markersize=5)    
        ax2.plot([minimum(t_coseismic), maximum(t_coseismic)], [3.464, 3.464], "k:", label="Vs(km/s)")  
        ax2.plot([minimum(t_coseismic), maximum(t_coseismic)], [5.996, 5.996], "k:", label="Vp(km/s)")    
        ax2.get_xaxis().set_tick_params(color=col)
        ax2.tick_params(axis="x", labelcolor=col)   
        # ax2.set_ylabel("Slip rate at the middle of nucleation zone(m/s)")
        ax2.set_ylabel("Rupture speed(km/s)")
        ax2.legend(loc="lower left")
        ax2.set_ylim([0, 7])

        # plot slip rate profile
        # ax = fig.add_subplot(n-1, 1, i)
        ax3 = fig.add_subplot(212)
        # println(size(t[indx_last_int:indx_last_int + N]))
        # println(size(value))
        
        println("The maximum moment release rate(Nm/s) is:", maximum(moment_rate))        
        ax3.plot(t_coseismic[2:end], moment_rate[2:end], color="red", label="moment rate")   
        # ax3.plot(0, moment_rate[nn], "*", color="black", markersize=15, label="nucleation to propagation")      
        ax3.set_xlabel("t(s)")
        ax3.set_yscale("log")
        # ax3.set_ylim([0, 10])
        ax3.set_ylabel("Moment rate(Nm/s)")
        ax3.legend(loc="upper left")
        ax3.set_title("Moment before propagation: $(round(moment_before, digits=2)) and corresponding magnitude: $(round(log10(moment_before)/1.5-6.07, digits=2))")

        col="tab:blue"
        ax4 = ax3.twinx()
        ax4.plot(t_coseismic[2:end], moment_acceleration[2:end], "o",color=col, label="moment acceleration") 
        # ax4.plot(0, moment_acceleration[nn], "*", color="black", markersize=15, label="Onset of rupture propagation")
        ax4.plot([0, 0], [maximum(moment_acceleration), minimum(moment_acceleration)], color="black")
        ax4.get_xaxis().set_tick_params(color=col)
        ax4.tick_params(axis="x", labelcolor=col) 
        ax4.set_yscale("log")
        ax4.set_ylabel("Moment acceleartion(Nm/s^{2})")
        ax4.legend(loc="lower left")

    end
    # println("Location and Full length of all seismic events' nucleation zone(km):", NS_width)
    # show()
    figname = string(path, "moment_release_example.png")
    fig.savefig(figname, dpi = 300)
end

function Nucleation(sliprate, FltX, tStart, t, N, criteria, measure_threshold)
    n_before = 100             # n_before*10 time steps before seismic threshold 
    n = length(tStart)         # how many seimsic events
    n_plot = 3      # number of plots of nucleation process
    NS_width = zeros(n-1,4)
    nn = 0
    plot_params()
    fig = PyPlot.figure(figsize=(10, 30))
    for i = 1: n-1     # neglect the first earthquake event
        
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:] .<= tStart[i+1])[end]   
        indx_last_int::Int = floor(indx_last/output_freq)
        # indx_last_int::Int = floor(indx_last)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        indx = findall(abs.(FltX) .<= 30)[1]
        value = sliprate[indx:end,indx_last_int-n_before:indx_last_int+N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
        depth = FltX[indx:end]

        for j = 2:N     # time step after the seismic threshold 1e-3m/s
            if maximum(value[:,n_before+j]) >= criteria
                nn = j
                break
            end
        end

        # measure the width of nucleation zone
        indx_nucleation = findall(value[:,n_before+nn] .>= measure_threshold)       # using the second line to define the width of nucleation size
        #println(indx_nucleation)
        if  indx_nucleation != []
            new_depth = FltX[indx:end][indx_nucleation]
            downdip_depth = maximum(new_depth)
            updip_depth = minimum(new_depth)

            NS_width[i, 1] = (downdip_depth + updip_depth)/2       # middle point 
            NS_width[i, 2] = downdip_depth - updip_depth        # width of nucleation zone
            NS_width[i, 3] = updip_depth                   # updip
            NS_width[i, 4] = downdip_depth          # downdip
        else
            continue
        end

        plot_params()   
        if i in n-n_plot:n-1
            # plot slip rate profile
            ax = fig.add_subplot(n_plot, 1, i-n+n_plot+1)
            # println(size(t[indx_last_int:indx_last_int + N]))
            # println(size(value))

            ax.plot(value[:,1:5:n_before+nn], depth, color="red")        # plot every five steps
            # ax.plot(value[:,2], depth, color="red")        # only plot the slip rate over seismicthreshold
            ax.set_xscale("log")
            ax.set_ylim([10,20])    
            ax.set_ylabel("Depth(km)")
            ax.set_xlim([1e-6, 1e-0])
            ax.set_xlabel("Slip Velocity(m/s)")
            ax.invert_yaxis()

            title = string(NS_width[i,2]," km")
            ax.set_title(title)
        end
    end
    # println("Location and Full length of all seismic events' nucleation zone(km):", NS_width)
    # show()
    figname = string(path, "sliprate_time_nucleation_alone.png")
    fig.savefig(figname, dpi = 300)
    return NS_width
end

function Nucleation_example(sliprate, weakeningrate, FltX, tStart, t, N, criteria, measure_threshold)
    n_before =100             # 200 time steps before seismic threshold 
    n = length(tStart)         # how many seimsic events
    # n = 5
    NS_width = zeros(n-1,4)
    nn = 0             # timestep when the maximum sliprate is over 1e-1m/s
    plot_params()
    fig = PyPlot.figure(figsize=(20, 10))
    # for i = 1: n-1 

    for i = n-3     # plot the i_th normal earthquake: choose by yourself    if i=n-1, then plot the last one!!
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:].<= tStart[i+1])[end]   
        indx_last_int::Int = floor(indx_last/output_freq)
        # indx_last_int::Int = floor(indx_last)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        indx = findall(abs.(FltX) .<= 30)[1]
        value = sliprate[indx:end,indx_last_int-n_before:indx_last_int+N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
        value_1 = weakeningrate[indx:end,indx_last_int-n_before:indx_last_int+N] 
        depth = FltX[indx:end]
        # find the timestep when sliprate first exceeds 1e-1 m/s
        for j = 2:N
            if  maximum(value[:,n_before+j]) >= criteria
                nn = j
                break
            end
        end
        println("The number of the timestep when maximum sliprate is over 0.1m/s is:", nn)
        # measure the width of nucleation zone
        indx_nucleation = findall(value[:,n_before + nn] .>= measure_threshold)       # using the second line(n_before+2) to define the width of nucleation size
        #println(indx_nucleation)
        new_depth = FltX[indx:end][indx_nucleation]
        downdip_depth = maximum(new_depth)
        updip_depth = minimum(new_depth)

        NS_width[i, 1] = (downdip_depth + updip_depth)/2       # middle point 
        NS_width[i, 2] = downdip_depth - updip_depth        # width of nucleation zone
        NS_width[i, 3] = updip_depth                   # updip
        NS_width[i, 4] = downdip_depth          # downdip

        # plot slip rate profile
        # ax = fig.add_subplot(n-1, 1, i)
        ax = fig.add_subplot(1, 2, 1)
        # println(size(t[indx_last_int:indx_last_int + N]))
        # println(size(value))
        
        ax.plot(depth, value[:,1:10:n_before+nn],color="red", )        # plot every five steps
        ax.set_yscale("log")
        ax.set_xlim([8,22])    
        ax.set_xlabel("Depth(km)")
        ax.set_ylim([1e-15, 1e1])
        ax.set_ylabel("Slip Velocity(m/s)")

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(depth, value_1[:,1:10:n_before+nn],color="blue", )        # plot every five steps
        ax2.set_yscale("log")
        ax2.set_xlim([8,22])    
        ax2.set_xlabel("Depth(km)")
        ax2.set_ylim([1e-10, 1e6])
        ax2.set_ylabel("weakeningrate")

        # ax.plot(value[:,1:1:2+n_before], depth, color="red", )        # plot every five steps
        # # ax.plot(value[:,2], depth, color="red")        # only plot the slip rate over seismicthreshold
        # ax.set_xscale("log")
        # ax.set_ylim([10,20])    
        # ax.set_ylabel("Depth(km)")
        # ax.set_xlim([1e-6, 1e-2])
        # ax.set_xlabel("Slip Velocity(m/s)")
        # ax.invert_yaxis()

        title = string(NS_width[i,2]," km")
        ax.set_title(title)
    end
    # println("Location and Full length of all seismic events' nucleation zone(km):", NS_width)
    # show()
    figname = string(path, "sliprate_time_nucleation_example.png")
    fig.savefig(figname, dpi = 300)
end


# Plot apparent_friction: with only every 0.1 second output
# function apparent_friction(stress, index_start, index_end, delfsec, index_ds_start,
#      index_ds_end, depth, t, N, normal_stress)

#     plot_params()
#     fig = PyPlot.figure(figsize=(7.2, 6))
#     ax = fig.add_subplot(111)

#     # average stress drop in nucleation_depth
#     seismic_stress = stress[index_start[N]:index_end[N], depth]
#     average_shear_stress_temp = mean(seismic_stress, dims = 2)
#     println(size(average_shear_stress_temp))
#     # average coseismic slip at nucleation_depth
#     seismic_delfsec = delfsec[index_ds_start[N]:index_ds_end[N], depth]
#     for i = 1:size(seismic_delfsec)[1]
#         seismic_delfsec[i,:] = seismic_delfsec[i,:] .- seismic_delfsec[1,:]
#     end
#     average_seismic_delfsec = mean(seismic_delfsec, dims = 2)

#     # resamlple the average shear stress to 0.1 s
#     t_real = t[10*index_start[N]: 10*index_end[N]]    # real timestep for average_shear_stress
#     len::Int = floor(t_real[end]-t_real[1])/0.1
#     average_shear_stress = zeros(len)
#     println(len)
#     j = 0
#     for i = 1:size(average_shear_stress_temp)[1]
#         if (t_real[i] - t_real[1]) >= 0.1*j
#             j = j+1
#             average_shear_stress[j] = average_shear_stress_temp[i]
#         endƒ
#     end

#     apparent_friction_coefficient = average_shear_stress ./ normal_stress

#     # calculate the minimum length
#     a = length(apparent_friction_coefficient)
#     b = length(average_seismic_delfsec)
#     len = minimum([a, b])

#     # println(apparent_friction_coefficient)
#     ax.plot(average_seismic_delfsec[1:len],apparent_friction_coefficient[1:len], lw = 2.0)
#     ax.set_xlabel("Slip(m)")
#     ax.set_ylabel("Apparent friction")
#     ax.set_ylim([0.0, 1.0])
#     ax.set_xlim([0.0, 6])
    
#     show()
    
#     figname = string(path, "apparent_friction.png")
#     fig.savefig(figname, dpi = 300)
# end


# Plot apparent_friction
function apparent_friction_new(stress, index_start, index_end, delfsec_et, index_ds_start,
    index_ds_end, NS_width, normal_stress, N_events)

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 6))
    ax = fig.add_subplot(111)
    print(size(NS_width)[1]) 
    for i = 1: N_events-1          # how many normal seismic events

        depth = findall(abs(NS_width[i,3]) .<= abs.(FltX) .<= abs(NS_width[i,4]))
        # average stress drop in nucleation_depth
        seismic_stress = stress[index_start[i+1]:index_end[i+1], depth]
        average_shear_stress = mean(seismic_stress, dims = 2)
        #println(size(average_shear_stress))
        #println(size(delfsec_et))
        # average coseismic slip at nucleation_depth
        seismic_delfsec = delfsec_et[index_ds_start[i+1]:index_ds_end[i+1], depth]   # the output frequency of coseismic slip is every 10 seconds
        #println(size(seismic_delfsec))
        # seismic_delfsec_temp = zeros(size(seismic_delfsec))

        # for j = 1:size(seismic_delfsec)[1]
        #     seismic_delfsec_temp[j,:] = seismic_delfsec[j,:] - seismic_delfsec[1,:]
        # end

        average_seismic_delfsec = mean(seismic_delfsec, dims = 2)
        #println("maximum value of displacement:", maximum(average_seismic_delfsec))
        apparent_friction_coefficient = average_shear_stress ./ normal_stress

        # calculate the minimum length
        a = length(apparent_friction_coefficient)
        #print(a)
        b = length(average_seismic_delfsec)
        #print(b)
        len = minimum([a, b])

        # println(apparent_friction_coefficient)

        ax.plot(average_seismic_delfsec[1:len], apparent_friction_coefficient[1:len], lw = 2.0,  label = i)

    end
    ax.set_xlabel("Slip(m)")
    ax.set_ylabel("Apparent friction")
    ax.set_ylim([0.4, 0.7])
    # ax.set_xlim([0.0, 6])
    ax.legend(loc="upper right") 
    
    show()
    
    figname = string(path, "apparent_friction.png")
    fig.savefig(figname, dpi = 300)
end

# Plot apparent_friction
function apparent_friction_new_prapogation(stress, index_start, index_end, delfsec_et, index_ds_start,
    index_ds_end, NS_width, normal_stress, N_events)

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 6))
    ax = fig.add_subplot(111)
    print(size(NS_width)[1]) 
    for i = 1: N_events-1      

        d = NS_width[i, 3] - 2    # 2 km shallower than the updip of the nucleation zone
        depth = findall(abs(d) - 0.5 .<= abs.(FltX) .<= abs(d) + 0.5)
        # average stress drop in nucleation_depth
        seismic_stress = stress[index_start[i+1]:index_end[i+1], depth]
        average_shear_stress = mean(seismic_stress, dims = 2)
        #println(size(average_shear_stress))
        #println(size(delfsec_et))
        # average coseismic slip at nucleation_depth
        seismic_delfsec = delfsec_et[index_ds_start[i+1]:index_ds_end[i+1], depth]   # the output frequency of coseismic slip is every 10 seconds
        #println(size(seismic_delfsec))
        # seismic_delfsec_temp = zeros(size(seismic_delfsec))

        # for j = 1:size(seismic_delfsec)[1]
        #     seismic_delfsec_temp[j,:] = seismic_delfsec[j,:] - seismic_delfsec[1,:]
        # end

        average_seismic_delfsec = mean(seismic_delfsec, dims = 2)
        #println("maximum value of displacement:", maximum(average_seismic_delfsec))
        apparent_friction_coefficient = average_shear_stress ./ normal_stress

        # calculate the minimum length
        a = length(apparent_friction_coefficient)
        #print(a)
        b = length(average_seismic_delfsec)
        #print(b)
        len = minimum([a, b])

        # println(apparent_friction_coefficient)

        ax.plot(average_seismic_delfsec[1:len], apparent_friction_coefficient[1:len], lw = 2.0,  label = i)

    end
    ax.set_xlabel("Slip(m)")
    ax.set_ylabel("Apparent friction")
    ax.set_ylim([0.4, 1.0])
    # ax.set_xlim([0.0, 6])
    ax.legend(loc="upper right") 
    
    show()
    
    figname = string(path, "apparent_friction_propagation_point.png")
    fig.savefig(figname, dpi = 300)
end




# Plot Vfmax
function velocity_dependence(b_value, Vfmax, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t./yr2sec, Vfmax, lw = 2.0, label="Max. Slip rate")
    lab1 = "Max. slip rate"
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    #ax.set_xlim([0, 600])
    ax.set_ylim([1e-10, 1e2])
    
    col="tab:red"
    ax2 = ax.twinx()
    
    ax2.plot(t./yr2sec, b_value, c=col, lw=2.0, label="b")
    lab2 = "b"
    ax.set_xlabel("Time (years)")
    ax2.set_ylabel("b of R&S friction law")
    ax2.set_ylim([0.015, 0.030])
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)

    #  ax.legend([lab1, lab2], loc=0)
    show()
    
    figname = string(path, "velocity_dependence_b.png")
    fig.savefig(figname, dpi = 300)
end



# Plot the stressdrops after each earthquake
function stressdrop_1(taubefore, tauafter, FltX)
    # N: number of subplots
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45));

    ax = fig.add_subplot(111)
      ax.plot(taubefore, FltX, lw = 2.0, color="tab:orange", 
              label="Shear stress before the earthquake", alpha=1.0);
      ax.plot(tauafter, FltX, lw = 2.0, color="tab:blue", 
              label="Shear stress after the earthquake", alpha=1.0);
      ax.set_xlabel("Stress drop (MPa)");
      ax.set_ylabel("Depth (km)");
      ax.set_ylim([0,5]);
      ax.set_xlim([5.5,6.5]);
      ax.invert_yaxis();
      plt.legend();
    
      show()
      figname = string(path, "shear_stress_1.png")
      fig.savefig(figname, dpi = 300)
    #  end
end



function cumSlipPlot_no_hypocenter(delfsec, delfyr, FltX, N)

    indx_1 = findall(abs.(FltX) .<= 25)[1]
    indx_2 = findall(abs.(FltX) .>= 5)[end]

    print(indx_1)
    print(indx_2)

    delfsec2 = transpose(delfsec[:, indx_1:indx_2])
    delfyr2 = transpose(delfyr)

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))

    ax = fig.add_subplot(111)
    plt.rc("font",size=12)
    ax.plot(delfyr2, FltX, color="royalblue", lw=1.0)
    ax.plot(delfsec2, FltX[indx_1:indx_2], color="chocolate", lw=1.0)
    # ax.plot(d_hypo, hypo./1000, "*", color="saddlebrown", markersize=25)
    ax.set_xlabel("Cumulative Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,30])

    L = (N+100) * 365 * 24 * 60 * 60 * 1e-9
    ax.set_xlim([0,L])
    #ax.set_xlim([0,9.0])
    ax.invert_yaxis()
    # show()
    
    figname = string(path, "cumulative_slip_no_hypocenter.png")
    fig.savefig(figname, dpi = 300)

end


# not yet working
function Nucleation_animation(sliprate, FltX, tStart, t, N, n)
    NS_width = zeros(n,1)
    plot_params()
    fig = PyPlot.figure(figsize=(5, 30))
for i = length(tStart)-n+1: length(tStart)
    #println("Time of the last seismic event(s):",tStart[end])
    indx_last = findall(t[:] .<= tStart[i])[end]   # last event!
    indx_last_int::Int = floor(indx_last/output_freq)
    #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

    indx = findall(abs.(FltX) .<= 20)[1]
    value = sliprate[indx:end,indx_last_int:indx_last_int+N]       # depth, timestep
    
    # NS_t_indx = findall(t .<= (tStart[i]+5))[end]  # what time(+15s) to measure the nucleation size
    # NS_t_indx_int::Int = floor(NS_t_indx/10)    # inde for output slip rate data
    # println("When to measure the nucleation size(s)", t[NS_t_indx]-t[indx_last])
    # println(NS_t_indx_int - indx_last_int)          # ~1 timestep
    
    # depth = FltX[indx:end]

    ax = fig.add_subplot(n,1, i+n-length(tStart))
    println(size(t[indx_last_int:indx_last_int + N]))
    println(size(value))
    
    # plot animation of nucleation process
    anim = @animate for i =1:10
        ax.plot(value[:,i], FltX[indx:end], color="blue" )
    end
    ax.set_xscale("log")
    ax.set_ylim([0,20])    
    ax.set_ylabel("Depth(km)")
    ax.set_xlim([1e-5, 1e1])
    ax.set_xlabel("Slip Velocity(m/s)")
    ax.invert_yaxis()
    #   cbar.set_ticks(cbar.get_ticks()[1:2:end])
    figname = string(path, "sliprate_time_",i,".gif")
    gif(anim, figname, fps=15)
end
#println("Full length of all seismic events' nucleation zone(km):", NS_width)
#show()

#fig.savefig(figname, dpi = 300)
end


# sliprate versus time for last 3 events
function eqCyclePlot_last_1(sliprate, FltX, tStart, t, N, n)
        NS_width = zeros(n,1)
        plot_params()
        fig = PyPlot.figure(figsize=(7.2, 10))
    for i = length(tStart)-n+1: length(tStart)
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:] .<= tStart[i])[end]   # last event!
        indx_last_int::Int = floor(indx_last/output_freq)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)
    
        indx = findall(abs.(FltX) .<= 20)[1]
        value = sliprate[indx:end,indx_last_int:indx_last_int+N]       # depth, timestep
        
        # NS_t_indx = findall(t .<= (tStart[i]+5))[end]  # what time(+15s) to measure the nucleation size
        # NS_t_indx_int::Int = floor(NS_t_indx/10)    # inde for output slip rate data
        # println("When to measure the nucleation size(s)", t[NS_t_indx]-t[indx_last])
        # println(NS_t_indx_int - indx_last_int)          # ~1 timestep
        
        # it is hard to record seismic events using the sliprate data
        # NS = findall(abs.(value[:, 1]).>= 1e-3) 
        # NS_indx_downdip = minimum(NS)
        # NS_indx_updip = maximum(NS)
        # updip_depth = FltX[NS_indx_updip]
        # downdip_depth = FltX[NS_indx_downdip]
        # NS_width[i-length(tStart)+n] = downdip_depth - updip_depth
        

        # depth = FltX[indx:end]
    
        # ax = fig.add_subplot(111)
        ax = fig.add_subplot(n,1, i+n-length(tStart))
        
        c = ax.imshow(value, cmap="magma", aspect="auto",
        #c = ax.imshow(value, cmap="inferno", aspect="auto",
                     #vmin=0.0, vmax=2.0,
                      norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e0),
                      interpolation="bicubic",
                      extent=[0, t[indx_last_int*10 + 10*N]- t[indx_last_int*10] , 0,20])
        
        # for stress
        #  c = ax.imshow(value, cmap="inferno", aspect="auto",
                      #  vmin=22.5, vmax=40,
                      #  interpolation="bicubic",
                      #  extent=[0,length(seismic_slipvel[1,:]), 0,16])
        ax.set_xlim([0,150])    
        ax.set_xlabel("Time(s)")
        ax.set_ylim([0, 20])
        ax.set_ylabel("Depth (km)")
        ax.invert_yaxis()
        cbar = fig.colorbar(c, label = "Slip rate(m/s)")
        #   cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    end
    show()
    figname = string(path, "sliprate_time_1.png")
    # figname = string(path, "sliprate_time_",i,".png")
    fig.savefig(figname, dpi = 300)
end

# sliprate versus time for last 3 events
function eqCyclePlot_last_2(sliprate, FltX, tStart, t,N,n)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 10))

    for i = length(tStart)-n+1: length(tStart)
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:] .<= tStart[i])[end]   # last event!
        indx_last_int::Int = floor(indx_last/output_freq)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)


        indx = findall(abs.(FltX) .<= 20)[1]
        value = sliprate[indx:end,indx_last_int:indx_last_int+N]       # depth, timestep
        
        # depth = FltX[indx:end]

        
        # ax = fig.add_subplot(111)
        ax = fig.add_subplot(n,1, i+n-length(tStart))
        
        #c = ax.imshow(value, cmap="magma", aspect="auto",
        c = ax.imshow(value, cmap="inferno", aspect="auto",
                    vmin=0.0, vmax=2.0,
                    #norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e0),
                    interpolation="bicubic",
                    extent=[0, t[indx_last_int*10 + 10*N]- t[indx_last_int*10] , 0,20])
        
        # for stress
        #  c = ax.imshow(value, cmap="inferno", aspect="auto",
                    #  vmin=22.5, vmax=40,
                    #  interpolation="bicubic",
                    #  extent=[0,length(seismic_slipvel[1,:]), 0,16])
            ax.set_xlim([0,150])    
            ax.set_xlabel("Time(s)")
            ax.set_ylim([0,20])
            ax.set_ylabel("Depth (km)")

        ax.invert_yaxis()
        cbar = fig.colorbar(c, label = "Slip rate(m/s)")
        #   cbar.set_ticks(cbar.get_ticks()[1:2:end])

    end
show()
figname = string(path, "sliprate_time_2.png")
# figname = string(path, "sliprate_time_",i,".png")
fig.savefig(figname, dpi = 300)
end




function velocity_dependence_b(x1 ,x2, y1, y2)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)            
    
    V::Vector{Float64} = range(1e-6, 1e-2, length=10000)
    b::Vector{Float64} = range(y1, y1, length=10000)        # y1 must be number
    b_linear::Vector{Float64} = range(y1, y1, length=10000)        # y1 must be number
    a::Vector{Float64} = range(0.015, 0.015, length=10000)        # y1 must be number

    index_1 = findall(x1.<= V .<= x2)
    index_2 = findall(x2.< V )
    b_linear[index_1] .= y1 .+ (V[index_1] .- x1)* (y2-y1) /(x2 - x1)
    # b[index_1] .= y1 .+ (log10.(V[index_1]) .- log10(x1)) * (y2-y1) / (log10(x2) - log10(x1))
    b[index_1] .= 0.019 .+ (0.025-0.019) .*  (1 .+sin.((V[index_1] .- 1e-5)./(1e-3 - 1e-5).*pi.-pi/2))./2
    b[index_2] .= y2
    ax.plot(V, b , "r", label="Evolution effect in R&S friction")
    ax.plot(V, b_linear , "r", linestyle=":", label="Evolution effect in R&S friction(linear)")
    ax.plot(V, a , "b", label="Direct effect in R&S friction")
    #ax.plot([1e-3, 1e-3],[0.01, 0.03] , "k", linestyle=":", label="Seismic threshold")
    ax.set_xlabel("Max. Slip rate (m/s)")
    ax.set_ylabel("b")
    ax.set_ylim([0.01, 0.03])
    ax.set_xlim([1e-5, 1e-3])
    # ax.set_xscale("log")
    ax.legend(loc="upper left") 
    show()
    figname = string(path, "velocity_dependence_b_seismogenic_zone.png")
    fig.savefig(figname, dpi = 300)

end


# depth vs event number
function hypo_Mw_stressdrop(hypo, Mw, del_sigma, delfafter, FltX)
  
    plot_params()
    fig = PyPlot.figure(figsize=(20,7))
    ax = fig.add_subplot(141)
    xaxis = hypo[Mw .>2.8]/1000    # choose seismic event larger than Mw2.8 

    ax.hist(xaxis[2:end], bins=collect(-0.5:16.5), orientation="horizontal", rwidth=0.9)
    ax.set_ylabel("Depth(km)")
    ax.set_xlabel("Event Number")
    ax.set_ylim([0,16])
    ax.invert_yaxis()

    ax = fig.add_subplot(142)
    xaxis = Mw[Mw .>2.8]   # choose seismic event larger than Mw2.8 

    ax.hist(xaxis[2:end], bins=collect(2.25:0.5:9.25), orientation="horizontal", rwidth=0.9)
    ax.set_ylabel("Moment Magnitude")
    ax.set_xlabel("Event Number")
    ax.set_ylim([2.5,9])
    ax.invert_yaxis()

    ax = fig.add_subplot(143)
    xaxis = del_sigma[Mw .>2.8]   # choose seismic event larger than Mw2.8 

    ax.hist(xaxis[2:end], bins=collect(0.05:0.2:5.05), orientation="horizontal", rwidth=0.9)
    ax.set_ylabel("Stress drop(MPa)")
    ax.set_xlabel("Event Number")
    # ax.set_ylim([1.5,3.5])
    ax.invert_yaxis()

    ax = fig.add_subplot(144)
    iter = length(delfafter[:,1])    # number of events
    xaxis = delfafter[Mw .>2.8,:]    # choose seismic event larger than Mw2.8 
    for i = 2:iter               # from the second earthquake
        ax.plot(xaxis[i,:], FltX)
        ax.set_xlabel("Coseismic Slip (m)")
        ax.set_ylabel("Depth (km)")
        ax.set_ylim([0,16])
        ax.set_xlim([0,maximum(delfafter)])
        #ax.set_xlim([0,9.0])
        ax.invert_yaxis()
    end
    
    show()
    figname = string(path, "hypo_Mw_stressdrop.png")
    fig.savefig(figname, dpi = 300)
end

# Plot shear stress comparison between full rupture and partial rupture in Prithvi(2021)
function shear_stress_comp(shear1b, shear1a, shear2b, shear2a, FltX1, FltX2)
    plot_params()
   
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    ax.plot(shear1b, FltX1, lw = 2.0, color="tab:blue",ls=:dashed, 
            label="Immature Fault Zone: before")
    ax.plot(shear1a, FltX1, lw = 2.0, color="tab:blue", label="Immature Fault Zone: after")
    ax.plot(shear2b, FltX2, lw = 2.0, color="tab:orange", ls=:dashed, 
            label="Mature Fault Zone: before")
    ax.plot(shear2a, FltX2, lw = 2.0, color="tab:orange", label="Mature Fault Zone: after")
    ax.set_xlabel("Shear stress (MPa)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,20])
    ax.invert_yaxis()
    plt.legend()
    show()
    
    figname = string(path, "Shear_Stress_002.png");
    fig.savefig(figname, dpi = 300);

end


# Plot rupture_length vs event time
# the effect of healing time??
function stem_plot(rl1, rl2, rl3, rl4)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.stem(-rl1./1e3, linefmt="C0-", markerfmt="C0o", basefmt=:None, label="Healing time = 10 yr")
    ax.stem(-rl2./1e3, linefmt="C1-", markerfmt="C1o", basefmt=:None, label="Healing time = 12 yr")
    ax.stem(-rl3./1e3, linefmt="C2-", markerfmt="C2o", basefmt=:None, label="Healing time = 15 yr")
    #  ax.stem(-rl4./1e3, linefmt="C3-", markerfmt="C3o", basefmt=:None, label="Healing time = 20 years")
    ax.set_xlabel("Event number")
    #  ax.set_xlabel("Event time (years)")
    ax.set_ylabel("Rupture length (km)")
    #  ax.set_xlim([10,300])
    plt.legend()
    show()
    
    figname = string(path, "stem_plot2.png")
    fig.savefig(figname, dpi = 300)

end

# Plot slip vs event number
function slipPlot(delfafter2, rupture_len, FltX, Mw, tStart)
    plot_params()
    fig, ax = PyPlot.subplots(nrows=1, ncols=4, sharex="all", sharey="all", figsize=(9.2, 5.00))

    xaxis = tStart[Mw .>2.8]    # choose seismic event larger than Mw2.8 
    println(xaxis)     # time for events
    delfafter = delfafter2[:,Mw .> 2.8]   # coseismic slip on fault for different events

    Mw2 = Mw[Mw .> 2.8]

    # Normalize colorbar
    norm = matplotlib.colors.Normalize(vmin = minimum(Mw2), 
                                       vmax=maximum(Mw2)) 
    colors = matplotlib.cm.inferno_r(norm(Mw2))

    #ax[1].barh([0,0.2,0.4,0.6,0.8,1], delfafter[end-1,:], height=6,        # coseismic slip at 60 m 
              #color=colors, align="center"); 
    ax[1].barh(xaxis, delfafter[end-1,:], height=6,        # coseismic slip at 60 m 
                color=colors, align="center");        # how to understand the function of xaxis
    # ax[1].set_ylabel("Time (yr)")
    # ax[1].invert_yaxis()
    # ax[1].set_title("At 60 m depth")

    # trench_depth1 = findall(abs.(FltX) .< 4.0e3)[1]
    # trench_depth2 = findall(abs.(FltX) .< 6.0e3)[1]
    # trench_depth3 = findall(abs.(FltX) .< 8.0e3)[1]
    
    # ax[2].barh(xaxis, delfafter[trench_depth1,:], height=6, 
    #           color=colors, align="center"); 
    # ax[2].invert_yaxis()
    # ax[2].set_title("At 4 km depth")
    
    # ax[3].barh(xaxis, delfafter[trench_depth2,:], height=6, 
    #           color=colors, align="center"); 
    # ax[3].invert_yaxis()
    # ax[3].set_title("At 6 km depth")
    
    # ax[4].barh(xaxis, delfafter[trench_depth3,:], height=6, 
    #           color=colors, align="center"); 
    # ax[4].invert_yaxis()
    # ax[4].set_title("At 8 km depth")
    
    # sm = matplotlib.cm.ScalarMappable(norm=norm, cmap="inferno_r")
    # sm.set_array([])
    # fig.colorbar(sm, shrink=0.9, label="Mw")
    # plt.xlabel("Coseismic Slip (m)")
    # plt.tight_layout()
    show()
    
    figname = string(path, "coseismic_slip.png")
    fig.savefig(figname, dpi = 300)
end


# Plot Vsurface
function VsurfPlot(Vsurf10, Vsurf12, Vsurf15, t10, t12, t15, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(8.2, 6.00))
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    ax.plot(t10./yr2sec, Vsurf10, lw = 2.0, label="Healing time = 10 yr")
    ax2.plot(t12./yr2sec, Vsurf12, lw = 2.0, label="Healing time = 12 yr")
    ax3.plot(t15./yr2sec, Vsurf15, lw = 2.0, label="Healing time = 15 yr")
    ax3.set_xlabel("Time (years)")
    ax2.set_ylabel("Surface. Slip rate (m/s)")
    ax.set_yscale("log")
    ax2.set_yscale("log")
    ax3.set_yscale("log")
    #  ax.set_xlim([230,400])
    #  plt.legend()
    plt.tight_layout()
    show()
    
    figname = string(path, "Vsurf02.png")
    fig.savefig(figname, dpi = 300)
end

# Plot alpha
function alphaaPlot(alphaa, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)

    ax.plot(t./yr2sec, alphaa.*100, lw = 2)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Shear Modulus Contrast (%)")
    #  ax.set_xlim([230,400])
    show()


    figname = string(path, "alpha_01.png")
    fig.savefig(figname, dpi = 300)
end

# Compare alpha
function alphaComp(a1, t1, a2, t2, a3, t3, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)

    ax.plot(t1./yr2sec, a1.*100, lw = 2, label="10 yr")
    ax.plot(t2./yr2sec, a2.*100, lw = 2, label="12 yr")
    ax.plot(t3./yr2sec, a3.*100, lw = 2, label="15 yr")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Shear Modulus Ratio (%)")
    #  ax.set_xlim([230,400])
    legend()
    show()


    figname = string(path, "alpha_comp.png")
    fig.savefig(figname, dpi = 300)
end

# Plot stressdrop comparison
function sd_comp(ds_im, tS_im, ds_m, tS_m)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.scatter(tS_im, ds_im, color="tab:blue", label="Immature")
    ax.scatter(tS_m, ds_m, color="tab:orange", label="Mature")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Stress drops (MPa)")
    plt.legend() 
    show()
    
    figname = string(path, "del_sigma_tStart.png")
    fig.savefig(figname, dpi = 300)
end