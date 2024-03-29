##############################
#  PLOTTING SCRIPTS
##############################

using Plots
using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
mpl = pyimport("matplotlib")
mpl.use("Agg")      

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
function icsPlot(a_b, Seff, tauo, FltX,Domain)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)            
    
    ax.plot(Seff, FltX, "b", label="Normal Stress")
    ax.plot(tauo, FltX, "orange", label="Shear Stress")
    ax.set_xlabel("Stress (MPa)")
    ax.set_ylabel("X (km)")
    ax.legend(loc="lower center") 
    
    col="tab:green"
    ax2 = ax.twiny()
    ax2.plot(a_b, FltX, "g",label="(a-b)")
    #println(FltX)
    # seismogenic_depth = findall(abs(Domain*Domain_X/4) .< abs.(FltX) .<= abs(Domain*Domain_X*3/4))   # note: unit of FltX here is km
    # a_b[seismogenic_depth] .= a_b[seismogenic_depth] .- 0.006
    #println(a_b)
    # ax2.plot(a_b, FltX, "r",label="coseismic (a-b)", linestyle=":")
    ax2.set_xlabel("Rate-state friction value (a-b)", color=col)
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)  
    ax2.set_xlim([-0.015,0.040])
    ax2.legend(loc="upper center") 

    ax.set_ylim([-Domain*Domain_X/2, Domain*Domain_X/2])
    ax.invert_yaxis()
             
    figname = string(path, "initial_condition.png")
    fig.savefig(figname, dpi = 300)
end


# Plot Vfmax
function ground_VPlot(ground_vel,Vfmax, out_seis_x, out_seis_y, N, t, criteria)
    

    plot_params()
    fig = PyPlot.figure(figsize=(15, 15))
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412)
    ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(414)

    N = 16000    # for expanding crack
    
    # t_seconds = (N+100) * 365 * 24 * 60 * 60 
    # indx_last = findall(t .<= t_seconds)[end]   # last event!
    indx::Int64 = findall(t[:] .<= tStart[5])[end] 
    # Start::Int = indx-N/2        # start time of the third event, Vmax>1e-3 m/s    however, it needs time for the wave or GPS to arrive at the station
    Start::Int = indx - N*4
    End::Int = indx + N*2
    # println("Epicenter distance(m):", out_seis_y)
    
    # println(minimum(Vfmax[indx:indx+N*2]))
    # println(maximum(Vfmax[indx:indx+N*2]))

    indx_nucleation::Int64 = indx + findall(Vfmax[indx:indx+N/2] .< criteria)[end]      # find the indx of Vf=1e-1 m/s    
    t_nucleation = t[indx_nucleation]
    println(indx_nucleation)

    indx_transition::Int64 = Start + findall(Vfmax[Start:indx+N/2] .< 1e-4)[end]      # find the indx of Vf=5e-3 m/s   scheme transition    
    # println(t[indx_nucleation])
    t_transition = t[indx_transition]
    println(indx_transition)

    # Start_1::Int = indx_transition - N/80        # -N/5 for fixed-length     -N/80 for expanding
    # End_1::Int = indx_nucleation + N/2        # +N for fixed length      +N/2   for expanding
    # t_win1 = t[Start_1]
    # t_win2 = t[End_1]

    for i = 5       # station number
        
        ax1.plot(t[Start: End], ground_vel[Start: End, i], lw = 2.0, label = string("Epicenter dist.(m): ", out_seis_y[i]))
        if i == 5
            ax1.plot([t_nucleation, t_nucleation],[-10, 10], "k", linestyle="-",linewidth=1, label="Time when on-fault sliprate is over 1e-1 m/s")
            ax1.plot([t_transition, t_transition],[-10, 10], "r", linestyle="--",linewidth=1, label="Time when on-fault sliprate is over 5e-3 m/s")
        else
            ax1.plot([t_nucleation, t_nucleation],[-10, 10], "k", linestyle="-",linewidth=1, )   
            ax1.plot([t_transition, t_transition],[-10, 10], "r", linestyle="--",linewidth=1, )
        end 

        ax2.plot(t[Start: End], ground_vel[Start: End, i], lw = 2.0,  label = string("Epicenter dist.(m): ", out_seis_y[i]))
        ax2.plot([t_nucleation, t_nucleation],[1e-15, 1e1], "k", linestyle="-",linewidth=1)
        ax2.plot([t_transition, t_transition],[1e-15, 1e1], "r", linestyle="--",linewidth=1)
        # plot the time window for ax3!
        # ax2.plot([t_win1, t_win1],[1e-15, 1e1], "r", linestyle="-",linewidth=1)
        # ax2.plot([t_win2, t_win2],[1e-15, 1e1], "r", linestyle="-",linewidth=1)

        # ax3.plot(t[Start_1: End_1], ground_vel[Start_1: End_1, i], lw = 1, linestyle="--",label = string("Epicenter dist.(m): ", out_seis_y[i]))
        # # ax3.plot([t_transition, t_transition],[-10, 10], "r", linestyle="--",linewidth=1)

        # # ax4.plot(t[Start+6000: End], ground_vel[Start+6000: End, i], lw = 2.0)
        # ax4.plot(t[Start_1: End_1], ground_vel[Start_1: End_1, i], lw = 1.0,linestyle="--",)
    end
    # ax.plot([0, indx_last],[1e-3, 1e-3] , "k", linestyle=":", label="Inertial threshold")
    # ax.plot([0, indx_last],[1e-1, 1e-1] , "k", linestyle="--", label="Seismic threshold")
    
    # ax2.legend(loc="upper left") 
    ax1.legend(loc="upper left") 
    ax1.set_xlabel("Time(s)")
    ax1.set_ylabel("Ground Velocity (m/s)")
    # ax2.set_yscale("log")
    ax1.set_ylim([-2, 5])         # here the velocity field could be negative!!

    ax2.set_xlabel("Time(s)")
    ax2.set_ylabel("Ground Velocity (m/s)")
    ax2.set_yscale("log")
    ax2.set_ylim([1e-10,1e1])        # m/s
    # ax.set_ylim([-2.5, 2.5])         # here the velocity field could be negative!!
    # show()

    # ax3.set_xlabel("Time(s)")
    # ax3.set_ylabel("Ground Velocity (m/s)")
    # # ax2.set_yscale("log")
    # ax3.set_ylim([-0.5e-3, 5e-3])
    # # ax3.set_xlim([t[Start_1], t[End_1]])
    # ax3.set_xlim([t[Start_1], t[Start_1]+10])      # only plot 10 seconds!

    # ax4.set_xlabel("Time(s)")
    # ax4.set_ylabel("Ground Velocity (m/s)")
    # ax4.set_yscale("log")
    # ax4.set_ylim([1e-6, 5e-3])  
    # # ax4.set_xlim([t_transition, t_nucleation])       #
    # ax4.set_xlim([t[Start_1], t[Start_1]+10])

    figname = string(path, "V_field.png")
    fig.savefig(figname, dpi = 300)

end

# Plot Vfmax
function VfmaxPlot(Vfmax, N, t)
    plot_params()
    fig = PyPlot.figure(figsize=(15, 5))
    ax = fig.add_subplot(111)

    t_seconds = (N+100) * 365 * 24 * 60 * 60 
    indx_last = findall(t .<= t_seconds)[end]   # last event!

    ax.plot(Vfmax[1:indx_last], lw = 2.0)
    ax.plot([0, indx_last],[1e-3, 1e-3] , "k", linestyle=":", label="Inertial threshold")
    ax.plot([0, indx_last],[1e-1, 1e-1] , "k", linestyle="--", label="Seismic threshold")
    ax.legend(loc="upper right") 
    ax.set_xlabel("Time steps")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.set_ylim([1e-10,1e2])
    # show()

    figname = string(path, "Vfmax.png")
    fig.savefig(figname, dpi = 300)
end

# Plot alpha and Vfmax on the same plot
function healing_analysis(Vf, alphaa, t, yr2sec, N)
    plot_params()
    fig = PyPlot.figure(figsize=(15, 5))
    ax = fig.add_subplot(111)
    
    ax.plot(t./yr2sec, Vf, lw = 2.0, label="Max. Slip rate")
    # ax.plot([0, maximum(t./yr2sec)],[1e-3, 1e-3] , "k", linestyle=":", label="Inertial threshold")
    ax.plot([0, maximum(t./yr2sec)],[1.3e-1, 1.3e-1] , "k", linestyle="--", label="Nucleation threshold")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.set_xlim([0, N])
    ax.set_ylim([1e-10, 1e2])
    ax.legend(loc="upper right") 

    # col="tab:red"
    # ax2 = ax.twinx()
    
    # ax2.plot(t./yr2sec, alphaa.*100, c=col, lw=2.0, label="Shear modulus ratio")
    # # lab2 = "Shear modulus ratio"
    # ax.set_xlabel("Time (years)")
    # ax2.set_ylabel("Shear Modulus (% of host rock)")
    # ax2.set_ylim([20, 110])
    # ax2.get_xaxis().set_tick_params(color=col)
    # ax2.tick_params(axis="x", labelcolor=col)

    #  ax.legend([lab1, lab2], loc=0)
    # show()
    
    figname = string(path, "healing_analysis.png")
    fig.savefig(figname, dpi = 300)
    for i in eachindex(Vf)
        if isnan(Vf[i]) == 1
            Vf[i] = 0
        end
    end
    n = findall(t .<= 10*yr2sec)[end]
    V_max = maximum(Vf[n: end])
    return V_max

 
end

function sliprate_analysis(Vf, alphaa, t, yr2sec, sliprate, FltX, N, Domain)

    plot_params()
    fig = PyPlot.figure(figsize=(8, 10))

    ax = fig.add_subplot(211)
    ax.plot(t./yr2sec, Vf, lw = 2.0)
    # ax.plot([0, maximum(t./yr2sec)],[1e-3, 1e-3] , "k", linestyle=":", label="Inertial threshold")
    ax.plot([0, maximum(t./yr2sec)],[1e-1, 1e-1] , "k", linestyle="--", label="Nucleation threshold")
    ax.set_xlabel("Year")
    ax.set_ylabel("Max. Slip Velocity (m/s)")
    ax.set_yscale("log")
    #ax.set_xlim([0, 600])
    ax.set_ylim([1e-10, 1e2])
    ax.set_xlim([N/2, N])
    # ax.set_xlim([0, N])
    ax.legend(loc="upper right") 
    ax.set_title("(a)", loc="left")

    t_seconds = N * 365 * 24 * 60 * 60 
    indx_last = findall(t .<= t_seconds)[end]   # last event!
    indx_last_int::Int = floor(indx_last/output_freq)

    indx = findall(abs.(FltX) .<= Domain_X*Domain)[1]
    value = sliprate[indx:end,1:indx_last_int]
    
    # print(findall(abs.(sliprate) .<= 1e-12))
    value[findall(abs.(value) .<= 1e-14)] .= 1e-14 
    
    # D = size(value)
    # y = FltX .* ones(D[2])'
    # println(size(y))
    # x = t[1:output_freq:end][1:D[2]]' .* ones(D[1]) # N
    # println(size(x))

    # depth = FltX[indx:end]

    t_seconds = N/2 * 365 * 24 * 60 * 60 
    indx_half = findall(t .<= t_seconds)[end]   # last event!
    indx_half_int::Int = floor(indx_half/output_freq)    # realistic timestep for sliprate data

    ax1 = fig.add_subplot(212)

    # c = ax1.pcolormesh(x, y, value, cmap="turbo",
    #             norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e-1),
    #             shading="gouraud",)

    c = ax1.imshow(value, cmap="turbo", aspect="auto",
                  norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e-1),
                  interpolation="bicubic",    # the interpolation method decide the final slip rate distrbution!!
                #   extent=[length(value[1,:])/2, length(value[1,:]), 0, Domain*Domain_X])
                  extent=[1,length(value[1,:]), 0, Domain*Domain_X])
    
    ax1.set_xlabel("Timestep")
    ax1.set_ylabel("X (km)")
    ax1.set_ylim([0, Domain*Domain_X])
    ax1.set_xlim([indx_half_int,length(value[1,:])])
    # ax1.set_xlim([length(value[1,:])/2,length(value[1,:])])

    # ax1.set_yticks([0, 2, 4, 6, 8, 10])
    # ax1.set_yticklabels([-5, -3, -1, 1, 3, 5])

    ax1.set_title("(b)", loc="left")
    ax1.invert_yaxis()
    
    fig.colorbar(c, ax=ax1, orientation = "horizontal", label="Slip Vlocity(m/s)")
    # fig.colorbar(c, ax=ax1, label="Slip rate(m/s)")
    
    # show()
    figname = string(path, "sliprate_analysis_1.png")
    fig.savefig(figname, dpi = 600)

end

function cumSlipPlot_no_hypocenter(delfsec, delfyr, FltX, N, Domain)

    indx_1 = findall(abs.(FltX) .<= 40)[1]
    indx_2 = findall(abs.(FltX) .>= 0)[end]

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
    ax.set_ylim([0,Domain*Domain_X])

    L = (N+100) * 365 * 24 * 60 * 60 * 1e-9
    ax.set_xlim([0,L])
    #ax.set_xlim([0,9.0])
    ax.invert_yaxis()
    # show()
    
    figname = string(path, "cumulative_slip_no_hypocenter.png")
    fig.savefig(figname, dpi = 300)

end

# Plot cumulative slip
function cumSlipPlot(delfsec, delfyr, FltX, hypo, d_hypo, N, Domain)
    
    indx_1 = findall(abs.(FltX) .<= Domain*Domain_X*7/8)[1]
    indx_2 = findall(abs.(FltX) .>= Domain*Domain_X/8)[end]

    print(indx_1)
    print(indx_2)

    delfsec2 = transpose(delfsec[:, indx_1:indx_2])
    delfyr2 = transpose(delfyr)

    plot_params()
    fig = PyPlot.figure(figsize=(20, 10))

    ax = fig.add_subplot(111)
    plt.rc("font",size=12)
    ax.plot(delfyr2, FltX, color="royalblue", lw=1.0)
    ax.plot(delfsec2, FltX[indx_1:indx_2], color="chocolate", lw=1.0)
    ax.plot(d_hypo, hypo./1000, "*", color="saddlebrown", markersize=25)
    ax.set_xlabel("Cumulative Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,Domain*Domain_X])

    L = (N+100) * 365 * 24 * 60 * 60 * 1e-9
    ax.set_xlim([0,L])
    #ax.set_xlim([0,9.0])
    ax.invert_yaxis()

    # show()
    
    figname = string(path, "cumulative_slip.png")
    fig.savefig(figname, dpi = 300)

end

# sliprate versus time plot
# N is the number of years to plot
function eqCyclePlot(sliprate, FltX, N, t, Domain)

    t_seconds = N * 365 * 24 * 60 * 60 
    indx_last = findall(t .<= t_seconds)[end]   # last event!
    indx_last_int::Int = floor(indx_last/output_freq)

    indx = findall(abs.(FltX) .<= Domain_X*Domain)[1]
    value = sliprate[indx:end,1:indx_last_int]
    
    # print(findall(abs.(sliprate) .<= 1e-12))
    value[findall(abs.(value) .<= 1e-14)] .= 1e-14 

    # depth = FltX[indx:end]

    plot_params()
    fig = PyPlot.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)

    c = ax.imshow(value, cmap="viridis", aspect="auto",
                  norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e-1),
                  interpolation="bicubic",    # the interpolation method decide the final slip rate distrbution!!
                  extent=[0,length(value[1,:]), 0, Domain*Domain_X])

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
    figname = string(path, "sliprate_time.png")
    fig.savefig(figname, dpi = 600)
    
end

function moment_release_example(sliprate, FltX, tStart, t, N, criteria, seismic_threshold, Domain)
    n_before = 150             # 200 time steps before seismic threshold 
    n = length(tStart)         # how many seimsic events
    # n = 5
    nn = 0             # timestep when the maximum sliprate is over 1e-1m/s
    plot_params()
    fig = PyPlot.figure(figsize=(10, 8))
    # for i = 1: n-1 
    
    for i = 5        # plot the i_th normal earthquake: choose by yourself    if i=n-1, then plot the last one!!
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:] .<= tStart[i+1])[end]   
        indx_last_int::Int = floor(indx_last/output_freq)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        t_coseismic = t[indx_last:output_freq:indx_last+N*output_freq]

        indx = findall(abs.(FltX) .<= Domain*Domain_X)[1]
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


        for k in 2:length(t_coseismic)     # time step after seismic threshold(1e-3 m/s)
                # measure the width of nucleation zone for each timestep
                # println("Max slip rate is", maximum(value[:, k]))
                indx_nucleation = findall(value[:, k] .>= seismic_threshold)       # using the second line(n_before+2) to define the width of nucleation size
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
        # ax.set_ylim([0, Domain*Domain_X/1.5])
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
        ax2.plot([minimum(t_coseismic), maximum(t_coseismic)], [5.996, 5.996], "k-", label="Vp(km/s)")    
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


function Nucleation(sliprate, weakeningrate, FltX, tStart, t, N, criteria, measure_threshold, Domain)
    n_before = 50             # n_before*10 time steps before seismic threshold 
    n = length(tStart)         # how many seimsic events
    n_plot = 3      # number of plots of nucleation process
    nn = 0
    mm = 0
    aa = 1
    NS_width = zeros(n-aa-2,4)
    min_Ω = zeros(n-aa-2)
    plot_params()
    fig = PyPlot.figure(figsize=(10, 30))
    for i = 1: n-aa-2     # neglect the last one event!  For some cases. this number could change to like n-10   
        
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:] .<= tStart[i+2])[end]   
        indx_last_int::Int = floor(indx_last/output_freq)
        # indx_last_int::Int = floor(indx_last)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)


        indx = findall(abs.(FltX) .<= Domain*Domain_X)[1]
        value = sliprate[indx:end,indx_last_int-n_before:indx_last_int+N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
        value_1 = exp.(weakeningrate[indx:end,indx_last_int-n_before:indx_last_int+N]).*value./1e-6 
        depth = FltX[indx:end]

        # for j = 2:N     # time step after the seismic threshold 1e-1m/s
        #     if maximum(value[:,n_before+j]) >= criteria
        #         nn = j
        #         break
        #     end
        # end

        nnn::Vector{Int64} = zeros(2)             # timestep when the maximum sliprate is over 1e-1m/s
        p = 1
        for j = 2:N
            if  (maximum(value[:,n_before+j]) >= criteria) && (maximum(value[:,n_before+j-1]) < criteria)
                nnn[p] = j
                p = p+1
                if p > 2
                    break     # only record 2 values
                end
            end
        end
        if nnn[2] != 0
            nn = nnn[1]        # only equals nnn[2] if a/b=0.9
        else
            nn = nnn[1]
        end
        println("Maximum sliprate is over 1e-1 m/s in these timesteps:", nnn)


        # find the index of maximum speed is over 1e-7m/s
        for j = 1: n_before     
            if maximum(value[:,j]) >= 1e-7
                mm = j
                break
            end
        end

        # measure the width of nucleation zone
        indx_nucleation = findall(value[:,n_before+nn] .>= measure_threshold)       # using the second line to define the width of nucleation size
        # println(indx_nucleation)
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

        # find the minimal Ω in the nucleation zone
        min_Ω[i] = minimum(value_1[indx_nucleation[:], mm:n_before+nn])

    end
    # # println("Location and Full length of all seismic events' nucleation zone(km):", NS_width)
    # # show()
    # figname = string(path, "sliprate_time_nucleation_alone.png")
    # fig.savefig(figname, dpi = 600)
    return NS_width, min_Ω       # must be in the end
end

function Nucleation_example(sliprate, weakeningrate, FltX, tStart, t, N, criteria, measure_threshold, Domain)
    n_before = 100             # 200 time steps before seismic threshold 
    n = length(tStart)         # how many seimsic events
    # n = 5
    NS_width = zeros(n-1,4)
    nn = 0             # timestep when the maximum sliprate is over 1e-1m/s
    mm = 0
    plot_params()
    fig = PyPlot.figure(figsize=(20, 30))
    # for i = 1: n-1 
    group = 3   # plot total 5 group of sliprate and weakenign rate
    for i = 1:group     # plot the i_th normal earthquake: choose by yourself    if i=n-1, then plot the last one!!
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:].<= tStart[i+1])[end]    # here i+1 is used
        indx_last_int::Int = floor(indx_last/output_freq)
        # indx_last_int::Int = floor(indx_last)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        indx = findall(abs.(FltX) .<= Domain*Domain_X)[1]
        value = sliprate[indx:end,indx_last_int-n_before:indx_last_int+N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
        value_1 = exp.(weakeningrate[indx:end,indx_last_int-n_before:indx_last_int+N]).*value./1e-6 
        depth = FltX[indx:end]

        # find the timestep when sliprate first exceeds 1e-1 m/s
        # for j = 2:N
        #     if  maximum(value[:,n_before+j]) >= criteria
        #         nn = j
        #         break
        #     end
        # end
        nnn::Vector{Int64} = zeros(2)             # timestep when the maximum sliprate is over 1e-1m/s
        p = 1
        for j = 2:N
            if  (maximum(value[indx_around_hypocenter,j]) >= criteria/(10^(n_e-k))) && (maximum(value[indx_around_hypocenter,j-1]) < criteria/(10^(n_e-k)))
                nnn[p] = j
                p = p+1
                if p > 2
                    break     # only record 2 values
                end
            end
        end
        if nnn[2] != 0
            nn = nnn[2]
        else
            nn = nnn[1]
        end
        println("Maximum sliprate is over 1e-1 m/s in these timesteps:", nnn)

        println("The number of the timestep when maximum sliprate is over 0.1m/s is:", nn)

        # find the index of maximum speed is over 1e-7m/s
        for j = 1: n_before     
            if maximum(value[:,j]) >= 1e-7
                mm = j
                break
            end
        end
            
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
        ax = fig.add_subplot(group, 2, 1+(i-1)*2)
        # println(size(t[indx_last_int:indx_last_int + N]))
        # println(size(value))
        
        ax.plot(depth, value[:,mm:5:n_before+nn],color="red", )        # plot every five steps
        ax.set_yscale("log")
        ax.set_xlim([Domain*Domain_X/8,Domain*Domain_X*7/8])   
        ax.set_xlabel("Depth(km)")
        ax.set_ylim([1e-15, 1e1])
        ax.set_ylabel("Slip Velocity(m/s)")

        ax2 = fig.add_subplot(group, 2, 2+(i-1)*2)
        ax2.plot(depth, value_1[:,mm:5:n_before+nn],color="blue", )        # plot every five steps
        ax2.plot([depth[1],depth[end]],[1, 1] , "k", linestyle=":", label="Ω=1")
        ax2.plot([depth[1],depth[end]],[2, 2] , "k", linestyle="--", label="Ω=2")
        ax2.plot([depth[1],depth[end]],[5, 5] , "k", linestyle="-", label="Ω=5")
        ax2.set_yscale("log")
        ax2.set_xlim([Domain*Domain_X/8,Domain*Domain_X*7/8])  
        ax2.set_xlabel("Depth(km)")
        ax2.set_ylim([1e-2, 1e10])
        ax2.set_ylabel("weakeningrate")
        ax2.legend(loc="upper right")

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

# only the characteristic events are plotted!!
function Nucleation_example_evolution(sliprate, weakeningrate, FltX, tStart, t, N, criteria, measure_threshold, Domain, TNS, hypo, characteristic_index)
    
    # time steps before seismic threshold 
    L = Domain*Domain_X/2
    println("The length of the fault(km):", L)
    println("The ratio of fault length to theretical value: ", L/TNS)

    if  0.5 <= L/TNS < 1  
        n_before = 400                    # 350, 450
    elseif 1 <= L/TNS < 2  
        n_before = 180        # remember to reduce this value when a/b is large!
        # n_before = 300
    # only works for a/b=0.9
    elseif 2 <= L/TNS < 4  
        n_before = 180    #  200 for case29,30,31,32,40     700 for case33,  900 for case34 and 35 in high_res_2
    elseif 4 <= L/TNS < 8 
        n_before = 130
    elseif 8 <= L/TNS < 16
        n_before = 80
    elseif 16 <= L/TNS <= 32
        n_before = 50
    else
        n_before = 30
    end

    n = length(tStart)         # how many seimsic events
    # n = 5
    # NS_width = zeros(n-1,4)

    nn = 0
    mm = 0
    plot_params()
    fig = PyPlot.figure(figsize=(20, 30))
    # for i = 1: n-1 

    group = length(characteristic_index)     # how many events will be displayed
    println("group: ", group)

    n_e = 7    # number of measurements of nulceation size
    NS_width = zeros(group, n_e, 4)
    EEP = zeros(group, n_e-1)
    EEPP = zeros(group, 1)

      # plot the i_th normal earthquake: choose by yourself    if i=n-1, then plot the last one!!
    for j = eachindex(characteristic_index)  
        #println("Time of the last seismic event(s):",tStart[end])
        i = characteristic_index[j] #
    # for i = 3:7
        indx_last = findall(t[:].<= tStart[i])[end]    # here i+2 is used: from the second normal event(or 3rd normal event)
        indx_last_int::Int = floor(indx_last/output_freq)
        # indx_last_int::Int = floor(indx_last)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        indx = findall(abs.(FltX) .<= Domain*Domain_X)[1]
        value = sliprate[indx:end,indx_last_int-n_before:indx_last_int+N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
        value_1 = exp.(weakeningrate[indx:end,indx_last_int-n_before:indx_last_int+N]).*value./1e-6    # weakening rate
        depth = FltX[indx:end]

        println("Mininum slip rate in selected data:", minimum(minimum(value[:,:])))
        println("Maxinum slip rate in selected data:", maximum(maximum(value[:,:])))

        indx_around_hypocenter = findall((hypo[i] .- TNS) .<= abs.(depth) .<= (hypo[i] .+ TNS))  # get the indx around nucleation location

        # println("Index around hypocenter: ", indx_around_hypocenter)
        # find the index of maximum speed is over 1e-8m/s
        for j = 1: n_before     
            if maximum(value[indx_around_hypocenter,j]) >= 1e-8
                mm = j
                break
            end
        end
        println("The number of the timestep over 1e-8 m/s is:", mm)

        for k = 1:n_e

            # find the timestep when sliprate first exceeds 1e-1 m/s, numbered from 1e-8 m/s
            # for j = 2:N
            #     if  maximum(value[indx_around_hypocenter,j]) >= criteria/(10^(n_e-k))
            #         nn = j
            #         break
            #     end
            # end

            # only for a/b = 0.9  
            nnn::Vector{Int64} = zeros(2)             # timestep when the maximum sliprate is over 1e-1m/s
            p = 1
            for j = 2:N
                if  (maximum(value[indx_around_hypocenter,j]) >= criteria/(10^(n_e-k))) && (maximum(value[indx_around_hypocenter,j-1]) < criteria/(10^(n_e-k)))
                    nnn[p] = j
                    p = p+1
                    if p > 2
                        break     # only record 2 values
                    end
                end
            end
            if nnn[2] != 0
                nn = nnn[1]    # here it should be nnn[2] for a/b=0.9
            else
                nn = nnn[1]
            end
            if k == n_e
                println("Maximum sliprate is over 1e-1 m/s in these timesteps:", nnn)
            end
            println("The number of the timestep from 1e-8 m/s when maximum sliprate is over ",criteria/(10^(n_e-k)) ,"m/s is:", nn - mm)

            # measure the width of nucleation zone

            indx_nucleation = findall(value[indx_around_hypocenter, nn] .>= criteria/(10^(n_e-k+1)))       # using the second line(n_before+2) to define the width of nucleation size
            # println(indx_nucleation)
            new_depth = new_depth = depth[indx_around_hypocenter][indx_nucleation]      # here the depth should be this!!
            downdip_depth = maximum(new_depth)
            updip_depth = minimum(new_depth)

            NS_width[j,k, 1] = (downdip_depth + updip_depth)/2       # middle point 
            NS_width[j,k, 2] = downdip_depth - updip_depth        # width of nucleation zone
            NS_width[j,k, 3] = updip_depth                   # updip
            NS_width[j,k, 4] = downdip_depth          # downdip

            
        end
        

        for k = 1:n_e-1
            EEP[j,k] = (log10(NS_width[j,k+1,2]*1e3) - log10(NS_width[j,k,2]*1e3))/(log10(criteria/(10^(n_e-k-1)))-log10(criteria/(10^(n_e-k))))
            # EEP[i,k] = (log10(NS_width[i,k+1,2]*1e3) - log10(NS_width[i,k,2]*1e3))/(log10(criteria/(10^(n_e-k-1)))-log10(criteria/(10^(n_e-k))))

        end

        EEPP[j] = sum(EEP[j,:])

        # plot slip rate profile
        # ax = fig.add_subplot(n-1, 1, i)
        ax = fig.add_subplot(group, 3, 1+(j-1)*3)
        # println(size(t[indx_last_int:indx_last_int + N]))
        # println(size(value))
        if 0.5 <= L/TNS < 1  
            n_inter = 30
        elseif 1 < L/TNS < 4  
            n_inter = 20 
        elseif 4< L/TNS < 16  
            n_inter = 10 
        else
            n_inter = 5
        end

        ax.plot(depth, value[:,mm:n_inter:nn],color="red", )        # plot every five steps
        ax.plot([depth[1],depth[end]],[3e-9, 3e-9] , "k", linestyle=":")
        ax.plot([depth[1],depth[end]],[2e-9, 2e-9] , "k", linestyle="--")
        ax.plot([depth[1],depth[end]],[1e-9, 1e-9] , "k", linestyle="-")
        ax.set_yscale("log")
        ax.set_xlim([Domain*Domain_X/8,Domain*Domain_X*7/8])   
        ax.set_xlabel("Depth(km)")
        ax.set_ylim([1e-15, 1e1])
        ax.set_ylabel("Slip Velocity(m/s)")
        
        ax2 = fig.add_subplot(group, 3, 2+(j-1)*3)
        ax2.plot(depth, value_1[:,mm:n_inter:nn],color="blue", )        # plot every five steps
        ax2.plot([depth[1],depth[end]],[1, 1] , "k", linestyle=":", label="Ω=1")
        ax2.plot([depth[1],depth[end]],[2, 2] , "k", linestyle="--", label="Ω=2")
        ax2.plot([depth[1],depth[end]],[5, 5] , "k", linestyle="-", label="Ω=5")
        ax2.set_yscale("log")
        ax2.set_xlim([Domain*Domain_X/8,Domain*Domain_X*7/8])  
        ax2.set_xlabel("Depth(km)")
        ax2.set_ylim([1e-2, 1e10])
        ax2.set_ylabel("weakeningrate")
        ax2.legend(loc="upper right")

        ax3 = fig.add_subplot(group, 3, 3+(j-1)*3)
        ax3.plot([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1], NS_width[j,:,2], color="limegreen", marker="o", label="Measured values")        # plot every five steps
        ax3.plot([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1], TNS.*ones(n_e,1), color="k", linestyle=":", label="Theoretical nucleation size")
        ax3.set_xlabel("Maximum velocity when make measurement(m/s)")
        ax3.set_xscale("log")
        ax3.set_ylabel("Nucleatin zone width(km)")
        ax3.legend(loc="upper right")

        # ax.plot(value[:,1:1:2+n_before], depth, color="red", )        # plot every five steps
        # # ax.plot(value[:,2], depth, color="red")        # only plot the slip rate over seismicthreshold
        # ax.set_xscale("log")
        # ax.set_ylim([10,20])    
        # ax.set_ylabel("Depth(km)")
        # ax.set_xlim([1e-6, 1e-2])
        # ax.set_xlabel("Slip Velocity(m/s)")
        # ax.invert_yaxis()

        title = string(NS_width[j,end,2]," km")
        ax.set_title(title)

        title = string("EEP=", EEPP[j])
        ax3.set_title(title)

    end
    
    # println("Location and Full length of all seismic events' nucleation zone(km):", NS_width)
    # show()
    
    figname = string(path, "sliprate_time_nucleation_example_evolution.png")
    fig.savefig(figname, dpi = 300)

    return mean(EEPP[:])
end

function Nucleation_example_no_weakening_rate(sliprate, FltX, tStart, t, N, criteria, measure_threshold)
    n_before = 100             # 200 time steps before seismic threshold 
    n = length(tStart)         # how many seimsic events
    # n = 5
    NS_width = zeros(n-1,4)
    nn = 0             # timestep when the maximum sliprate is over 1e-1m/s
    plot_params()
    fig = PyPlot.figure(figsize=(20, 10))
    # for i = 1: n-1 

    for i = n-5     # plot the i_th normal earthquake: choose by yourself    if i=n-1, then plot the last one!!
        #println("Time of the last seismic event(s):",tStart[end])
        indx_last = findall(t[:].<= tStart[i+1])[end]   
        indx_last_int::Int = floor(indx_last/output_freq)
        # indx_last_int::Int = floor(indx_last)
        #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        indx = findall(abs.(FltX) .<= 30)[1]
        value = sliprate[indx:end,indx_last_int-n_before:indx_last_int+N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
        # value_1 = exp.(weakeningrate[indx:end,indx_last_int-n_before:indx_last_int+N]).*value./1e-6 
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
        ax = fig.add_subplot(1, 1, 1)
        # println(size(t[indx_last_int:indx_last_int + N]))
        # println(size(value))
        
        ax.plot(depth, value[:,1:5:n_before+nn],color="red", )        # plot every five steps
        ax.set_yscale("log")
        ax.set_xlim([10,20])    
        ax.set_xlabel("Depth(km)")
        ax.set_ylim([1e-15, 1e1])
        ax.set_ylabel("Slip Velocity(m/s)")

        # ax2 = fig.add_subplot(1, 2, 2)
        # ax2.plot(depth, value_1[:,1:5:n_before+nn],color="blue", )        # plot every five steps
        # ax2.plot([depth[1],depth[end]],[1, 1] , "k", linestyle=":", label="Ω=1")
        # ax2.set_yscale("log")
        # ax2.set_xlim([8,22])    
        # ax2.set_xlabel("Depth(km)")
        # ax2.set_ylim([1e-2, 1e10])
        # ax2.set_ylabel("weakeningrate")

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
    figname = string(path, "sliprate_time_nucleation_example_no_weakening_rate.png")
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
#         end
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
function stressdrop_2(taubefore, tauafter, FltX, tStart, Domain)
    N = length(tStart)
    nn = 3
    plot_params()
    fig = PyPlot.figure(figsize=(10, 20));
    for i = 1:nn
    ax = fig.add_subplot(nn,1,i)
      ax.plot(taubefore[i,:], FltX, lw = 2.0, color="tab:orange", 
              label="Shear stress before the earthquake", alpha=1.0);
      ax.plot(tauafter[i,:], FltX, lw = 2.0, color="tab:blue", 
              label="Shear stress after the earthquake", alpha=1.0);
      ax.set_xlabel("Stress drop (MPa)");
      ax.set_ylabel("Depth (km)");
      ax.set_ylim([0,Domain*Domain_X]);
      ax.set_xlim([10,45]);
      ax.invert_yaxis();
      plt.legend();
    end
      show()
      figname = string(path, "shear_stress_following.png")
      fig.savefig(figname, dpi = 300)
    #  end
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
      ax.set_ylim([0,20]);
      ax.set_xlim([15,45]);
      ax.invert_yaxis();
      plt.legend();
    
      show()
      figname = string(path, "shear_stress_1.png")
      fig.savefig(figname, dpi = 300)
    #  end
end






# not yet working
function Nucleation_animation(sliprate, FltX, tStart, t, N, n)
    NS_width = zeros(n,1)
    plot_params()
    fig = PyPlot.figure(figsize=(5, 30))
for i = length(tStart)-n+1: length(tStart)
    #println("Time of the last seismic event(s):",tStart[end])
    indx_last = findall(t[:] .<= tStart[i])[end]   # last event!
    indx_last_int::Int = floor(indx_last/10)
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
function eqCyclePlot_last_1(sliprate, FltX, tStart, t, N, n, Domain, characteristic_index, TNS, vs2)
        plot_params()
        fig = PyPlot.figure(figsize=(12, 20))

    for j in eachindex(characteristic_index)   # how many selected characteristic events
        #println("Time of the last seismic event(s):",tStart[end])
        i = characteristic_index[j]      # real index of characteristic events
        indx_last = findall(t[:] .<= tStart[i])[end]  
        indx_last_int::Int = floor(indx_last/output_freq)
        println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        value = sliprate[:,indx_last_int:indx_last_int+N]       # depth, timestep
        
        # If the minimum sliprate is too small, the interpolation will fail, 
        # so we need to neglect the slip rate below 1e-12 as following
        # print(findall(abs.(sliprate) .<= 1e-12))
        value[findall(abs.(value) .<= 1e-14)] .= 1e-14 

        D = size(value)
        println("Maximum value of slip rate: ", maximum(value))
        println("Minimum value of slip rate: ", minimum(value))
        println(D)
        y = FltX .* ones(D[2])'
        println(size(y))
        x = t[indx_last_int*output_freq:10:indx_last_int*output_freq + output_freq*N]' .* ones(D[1]) # N
        println(size(x))

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
        ax = fig.add_subplot(n/2,2, j)


        # c = ax.imshow(value, cmap="jet", aspect="auto",
        #              #vmin=0.0, vmax=2.0,
        #               norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e1),
        #               interpolation="bicubic",
        #               extent=[0, t[indx_last_int*output_freq + output_freq*N]- t[indx_last_int*output_freq], 0, Domain*Domain_X])
        
        c = ax.pcolormesh(x, y, value, cmap="turbo",
        #c = ax.imshow(value, cmap="inferno", aspect="auto",
                    #vmin=0.0, vmax=2.0,
                    # interpolation="bicubic",
                    norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e1),
                    shading="gouraud",)
              

        # for stress
        #  c = ax.imshow(value, cmap="inferno", aspect="auto",
                      #  vmin=22.5, vmax=40,
                      #  interpolation="bicubic",
                      #  extent=[0,length(seismic_slipvel[1,:]), 0,16])
        # ax.set_xlim([0, 1e4])        # unit: seconds
        ax.set_xlabel("Time(s)")
        ax.set_ylim([0, Domain*Domain_X])
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
function eqCyclePlot_last_2(sliprate, FltX, tStart, t, N, n, Domain, characteristic_index, TNS, vs2)
    plot_params()
    fig = PyPlot.figure(figsize=(12, 20))

    for j in eachindex(characteristic_index)   # how many selected characteristic events
        #println("Time of the last seismic event(s):",tStart[end])
        i = characteristic_index[j]      # real index of characteristic events
        indx_last = findall(t[:] .<= tStart[i])[end]  
        indx_last_int::Int = floor(indx_last/output_freq)
        println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        value = sliprate[:,indx_last_int:indx_last_int+N]       # depth, timestep
        
        # If the minimum sliprate is too small, the interpolation will fail, 
        # so we need to neglect the slip rate below 1e-12 as following
        # print(findall(abs.(sliprate) .<= 1e-12))
        value[findall(abs.(value) .<= 1e-14)] .= 1e-14 

        D = size(value)
        println("Maximum value of slip rate: ", maximum(value))
        println("Minimum value of slip rate: ", minimum(value))
        println(D)
        y = FltX .* ones(D[2])'
        println(size(y))
        x = t[indx_last_int*output_freq:10:indx_last_int*output_freq + output_freq*N]' .* ones(D[1]) # N
        println(size(x))
        
        ax = fig.add_subplot(n/2,2, j)
        c = ax.pcolormesh(x, y, value, cmap="seismic",
            norm=matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e1),
            shading="gouraud",)
        ax.set_xlabel("Time(s)")
        ax.set_ylim([0, Domain*Domain_X])
        ax.set_ylabel("Depth (km)")
        ax.invert_yaxis()
        cbar = fig.colorbar(c, label = "Slip rate(m/s)")   
        
    end

figname = string(path, "sliprate_time_2.png")
# figname = string(path, "sliprate_time_",i,".png")
fig.savefig(figname, dpi = 300)
end

# sliprate versus time for last 3 events
function eqCyclePlot_last_3(sliprate, FltX, tStart, t, N, n, Domain, characteristic_index, TNS, vs2)
    plot_params()
    fig = PyPlot.figure(figsize=(12, 20))

    for j in eachindex(characteristic_index)   # how many selected characteristic events
        #println("Time of the last seismic event(s):",tStart[end])
        i = characteristic_index[j]      # real index of characteristic events
        indx_last = findall(t[:] .<= tStart[i])[end]  
        indx_last_int::Int = floor(indx_last/output_freq)
        println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)

        value = sliprate[:,indx_last_int:indx_last_int+N]       # depth, timestep
        
        # If the minimum sliprate is too small, the interpolation will fail, 
        # so we need to neglect the slip rate below 1e-12 as following
        # print(findall(abs.(sliprate) .<= 1e-12))
        value[findall(abs.(value) .<= 1e-14)] .= 1e-14 

        D = size(value)
        println("Maximum value of slip rate: ", maximum(value))
        println("Minimum value of slip rate: ", minimum(value))
        println(D)
        # y = FltX .* ones(D[2])'
        # println(size(y))
        # x = t[indx_last_int*output_freq:10:indx_last_int*output_freq + output_freq*N]' .* ones(D[1]) # N
        # println(size(x))
        
        ax = fig.add_subplot(n/2,2, j)

        c = ax.imshow(value, cmap="seismic", aspect="auto",
                norm=matplotlib.colors.LogNorm(vmin=1e-1, vmax=1e1),
                interpolation="bicubic",    # the interpolation method decide the final slip rate distrbution!!
                extent=[0,length(value[1,:]), 0, Domain*Domain_X])
        
        ax.set_xlabel("Time step")
        ax.set_ylim([0, Domain*Domain_X])
        ax.set_ylabel("Depth (km)")
        ax.invert_yaxis()
        cbar = fig.colorbar(c, label = "Slip rate(m/s)")   
        
    end

figname = string(path, "sliprate_time_3.png")
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
    b[index_1] .= 0.019 .+ (0.025 - 0.019) .* (1 .+sin.((V[index_1] .- 1e-5)./(1e-3 - 1e-5).*pi.-pi/2))./2
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
function hypo_Mw_stressdrop(hypo, Mw, del_sigma, delfafter, FltX, Domain)
  
    plot_params()
    fig = PyPlot.figure(figsize=(20,7))
    ax = fig.add_subplot(141)
    xaxis = hypo[Mw .>2.8]/1000    # choose seismic event larger than Mw2.8 

    ax.hist(xaxis[2:end], bins=collect(-0.5:16.5), orientation="horizontal", rwidth=0.9)
    ax.set_ylabel("Depth(km)")
    ax.set_xlabel("Event Number")
    ax.set_ylim([0,Domain*Domain_X])
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
    println("stress drops(MPa)", xaxis)

    ax.hist(xaxis[2:end], bins=collect(0.05:0.2:15.05), orientation="horizontal", rwidth=0.9)
    ax.set_ylabel("Stress drop(MPa)")
    ax.set_xlabel("Event Number")
    # ax.set_ylim([1.5,3.5])
    ax.invert_yaxis()

    ax = fig.add_subplot(144)
    xaxis = delfafter[Mw .>2.8,:]    # choose seismic event larger than Mw2.8 
    iter = length(xaxis[:,1])    # number of events
    println("Mw>2.8 events number:", iter)
    for i = 2:iter               # from the second earthquake
        ax.plot(xaxis[i,:], FltX)
        ax.set_xlabel("Coseismic Slip (m)")
        ax.set_ylabel("Depth (km)")
        ax.set_ylim([0,Domain*Domain_X])
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