using DelimitedFiles
using LinearAlgebra

# include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")

multiple = [4.0,5.0,6.0,7.0]
cos_reduction = [0.05, 0.06,0.07,0.08]
title = ["(a)", "(b)", "(c)", "(d)"]

# Global variables
yr2sec = 365*24*60*60
# comment this part if there is nothing in event_time temporarily
fig = PyPlot.figure(figsize=(6, 15));
fig.patch.set_facecolor("#fff8e5")
plt.rc("font",size=12)

for i = 1:4                   # normal stress
    # for j = 1:4              # cos_reduction 

        FILE = "20000.0_1500.0_8_0.85_$(cos_reduction[i])_$(multiple[i])"   # normal stress testing
        println(FILE)
        out_path = "$(@__DIR__)/data/immature_fully_healing/$(FILE)/"


                # Order of storage: Seff, tauo, FltX, cca, ccb, xLf
                params = readdlm(string(out_path, "params.out"), header=false)

                # Seff = params[1,:]
                # tauo = params[2,:]
                FltX = params[3,:]
                # cca = params[4,:]
                # ccb = params[5,:]
                # a_b = cca .- ccb
                # Lc = params[6,:]

        # # Read data
        # event_time = readdlm(string(out_path, "event_time.out"), header=false)
        # tStart = event_time[:,1]
        # #println(tStart)
        # tEnd = event_time[:,2]

        # hypo = event_time[:,3]
        # d_hypo = event_time[:,4]    # unit: m 
        # println(d_hypo) 
        # println(hypo)


        # displacement on fault line for different time 
        delfsec = readdlm(string(out_path, "delfsec.out"))
        # print(size(delfsec))
        delfyr = readdlm(string(out_path, "delfyr.out"))
        # print(size(delfyr))

        # culmulative slip
        #cumSlipPlot(delfsec[1:end,:], delfyr[1:end, :], FltX);
        indx = findall(abs.(FltX) .<= 15)[1]           
        delfsec2 = transpose(delfsec[1:5:end,indx:end])  # 0.5 seconds
        delfyr2 = transpose(delfyr[1:2:end,:])       # 4 year 

        delfsec = 0
        delfyr = 0

        ax = fig.add_subplot(4,1, i)
        plt.rc("font",size=12)
        ax.plot(delfyr2, FltX, color="royalblue", lw=1.0)
        ax.plot(delfsec2, FltX[indx:end], color="chocolate", lw=1.0)
        #ax.plot(d_hypo, hypo./1000 , "*", color="saddlebrown", markersize=20)
        ax.set_xlabel("Cumulative slip (m)")
        ax.set_ylabel("Depth (km)")
        ax.set_ylim([0,20])
        ax.set_xlim([0,maximum(delfyr2)])
        ax.invert_yaxis()
        ax.set_title(title[i], loc="left")
        # # # healing analysis: Vfmax and regidity ratio vs. time
        #healing_analysis(Vfmax, alphaa, t, yr2sec)
               
        delfsec2 = 0
        delfyr2 = 0
        
end
path = "$(@__DIR__)/plots/immature_fully_healing/"        
figname = string(path, "cumulative_slip_combined.png")
fig.savefig(figname, dpi = 600, bbox_inches="tight", pad_inches=0)
#show()



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