using DelimitedFiles
using LinearAlgebra

# include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")

Faultzone_depth::Vector{Int64} = [6000, 8000, 10000, 12000]
title = ["(a)", "(b)", "(c)", "(d)"]

output_freq = 10

# Global variables
yr2sec = 365*24*60*60
# comment this part if there is nothing in event_time temporarily
fig = PyPlot.figure(figsize=(25, 10));
# fig.patch.set_facecolor("#fff8e5")

mpl = pyimport("matplotlib")


for i = 1:2                  # normal stress
    # for j = 1:4              # cos_reduction 

        FILE = "1.0_10_1000_$(Faultzone_depth[i])_150_0.008_5"   # normal stress testing
        println(FILE)
        out_path = "$(@__DIR__)/data/benchmark_project/$(FILE)/"


                # Order of storage: Seff, tauo, FltX, cca, ccb, xLf
                params = readdlm(string(out_path, "params.out"), header=false)

                # Seff = params[1,:]
                # tauo = params[2,:]
                FltX = params[3,:]
                # cca = params[4,:]
                # ccb = params[5,:]
                # a_b = cca .- ccb
                # Lc = params[6,:]

        time_vel = readdlm(string(out_path, "time_velocity.out"), header=false)
        t = time_vel[:,1]             # all real timsteps

        # Read data
        event_time = readdlm(string(out_path, "event_time.out"), header=false)
        tStart = event_time[:,1]
        #println(tStart)
        tEnd = event_time[:,2]

        # hypo = event_time[:,3]
        d_hypo = event_time[:,4]    # unit: m 
        # println(d_hypo) 
        # println(hypo)

        # find the hypocenter
        sliprate = readdlm(string(out_path, "sliprate.out"), header=false)   # every 10 timesteps
        sliprate = sliprate'
        println("Dimension of sliprate:",size(sliprate))

        n = length(tStart)         # how many seimsic events
        hypo::Vector{Float64} = zeros(n)

        mm = 0
        N = 1000
    
          # plot the i_th normal earthquake: choose by yourself    if i=n-1, then plot the last one!!
        for i = 1:n

            indx_last = findall(t[:].<= tStart[i])[end]    # here i+2 is used: from the second normal event(or 3rd normal event)
            indx_last_int::Int = floor(indx_last/output_freq)
            # indx_last_int::Int = floor(indx_last)
            #println("Index of timestep in sliprate(output every 10) at the beginning of last seismic event:", indx_last_int)
    
            indx = findall(abs.(FltX) .<= 24)[1]
            value = sliprate[indx:end,indx_last_int:indx_last_int + N]       # depth, timestep  # only use the slip rate data after the seismic threshold!
            # value_1 = exp.(weakeningrate[indx:end,indx_last_int-n_before:indx_last_int+N]).*value./1e-6    # weakening rate
            depth = FltX[indx:end]
    
            # println("Index around hypocenter: ", indx_around_hypocenter)
            # find the index of maximum speed is over 1e-8m/s
            for j = 1: N     
                if maximum(value[:,j]) >= 5e-3
                    mm = j
                    break
                end
            end
            println("The number of the timestep over 5e-3 m/s is:", mm)

            indx_nucleation = findall(value[:, mm] .>= 1e-3)       # using the second line(n_before+2) to define the width of nucleation size
            println(indx_nucleation)
            new_depth = new_depth = depth[indx_nucleation]      # here the depth should be this!!
            downdip_depth = maximum(new_depth)
            updip_depth = minimum(new_depth)

            hypo[i] = (downdip_depth + updip_depth)/2       # middle point 
        end

        # displacement on fault line for different time 
        delfsec = readdlm(string(out_path, "delfsec.out"))
        # print(size(delfsec))
        delfyr = readdlm(string(out_path, "delfyr.out"))
        # print(size(delfyr))

        # culmulative slip
        #cumSlipPlot(delfsec[1:end,:], delfyr[1:end, :], FltX);
        indx = findall(abs.(FltX) .<= 24)[1]           
        delfsec2 = transpose(delfsec[1:10:end,indx:end])  # 1 seconds
        delfyr2 = transpose(delfyr[1:10:end,:])       # 20 year 

        delfsec = 0
        delfyr = 0

        ax = fig.add_subplot(2, 8, (4*(i-1)+1,4*(i-1)+3))
        ax.plot(delfyr2, FltX, color="royalblue", lw=1.0)
        ax.plot(delfsec2, FltX[indx:end], color="chocolate", lw=1.0)
        ax.plot(d_hypo, hypo , "*", color="saddlebrown", markersize=20)
        ax.set_xlabel("Cumulative slip (m)")
        ax.set_ylabel("Depth (km)")
        ax.set_ylim([0,24])
        ax.set_xlim([0,maximum(delfyr2)])
        ax.invert_yaxis()
        ax.set_title(string(Faultzone_depth[i]./1000, "km"), loc="center")
        # # # healing analysis: Vfmax and regidity ratio vs. time
        #healing_analysis(Vfmax, alphaa, t, yr2sec)

        ax2 = fig.add_subplot(2, 8, 4*i)
        ax2.hist(hypo, bins=collect(-0.5:24.5), density="True", orientation="horizontal", rwidth=0.9, )
        ax2.set_ylabel("Depth(km)")
        ax2.set_xlabel("Frequency(%)")
        ax2.set_ylim([0,24])
        ax2.invert_yaxis()
        ax2.xaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1))
        plt.tight_layout()   # adjust the spacing between subplots to minimize the overlaps
        
end
path = "$(@__DIR__)/plots/benchmark_project/"        
figname = string(path, "cumulative_slip_combined.png")
fig.savefig(figname, dpi = 600, bbox_inches="tight", pad_inches=0)
# fig.savefig(figname, dpi = 600)
#show()