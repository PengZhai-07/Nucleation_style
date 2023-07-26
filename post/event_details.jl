#################################
# MODULE FOR SOME CALCULATIONS
# FROM SIMULATION OUTPUT
#################################
using LinearAlgebra

#..........................................................
# Compute the moment magnitude:
#       Assumed the rupture area to be square; the rupture
#       dimension along depth is the same as the rupture
#       dimension perpendicular to the plane
#..........................................................
function moment_magnitude_new(mu, FltX, delfafter, stressdrops, delfsec, index_ds_start, index_ds_end, stress, index_start, index_end)
    # Final coseismic slip of each earthquake
    FltNglob = length(FltX)         # number of GLL nodes on fault

    iter = length(delfafter[1,:])    # number of events
    Mw = zeros(iter)
    seismic_moment = zeros(iter)
    rupture_len = zeros(iter)
    fault_slip = zeros(iter)
    E_T0 = zeros(iter)
    E_R = zeros(iter)
    ER_M0_ratio = zeros(iter)
    radiation_eff = zeros(iter)

    del_sigma = zeros(iter)
    
    # calculate the nearest grid size
    dx = diff(FltX).*1e3        # Fltx: real depth(km) of all nodes on the fault 

    for i = 1:iter-1    # for each recorded event(including the first artificial event)
        
        # slip threshold = 1% of maximum final slip
        slip_thres = 0.01*maximum(delfafter[:,i])

        # area = slip*(rupture dimension along depth)
        # zdim = rupture along z dimension = depth rupture dimension
        area = 0; zdim = 0; temp_sigma = 0; 
        temp_sigma = 0
        temp_E_T0 = 0
        temp_E_R = 0

        for j = 1:FltNglob      # total points on fault
            if  delfafter[j,i] >= slip_thres    # only the points with final slip > 1% of maximum slip are counted
                area = area + delfafter[j,i]*dx[j-1]
                zdim = zdim + dx[j-1]     

                # Avg. stress drops along rupture area
                temp_sigma = temp_sigma + stressdrops[j,i]*dx[j-1]

                # available energy
                temp_E_T0 = temp_E_T0 + stressdrops[j,i]*delfafter[j,i]*(-dx[j-1])    # stress(MPa) * slip(m)

                n = index_ds_end[i] - index_ds_start[i]
                # println("n")
                temp_E_R = temp_E_R + sum((stress[index_start[i]:index_start[i]+n-1,j].-stress[index_start[i]+n,j]).*(delfsec[index_ds_start[i]+1:index_ds_end[i],j].- delfsec[index_ds_start[i]:index_ds_end[i]-1,j]))*(-dx[j-1])
                # temp_E_R = temp_E_R + sum((stress[index_start[i]:index_end[i]-1,j].-stress[index_end[i],j]).*(delfsec[index_ds_start[i]+1:index_ds_end[i],j].- delfsec[index_ds_start[i]:index_ds_end[i]-1,j]))*(-dx[j-1])
            end
        end

        E_T0[i] = 0.5*temp_E_T0   # based on Huang et al.(2014) equation (2)
        # println("available energy(MPa*m^2):",E_T0[i])
        # println("fracture(rupture energy(MPa*m^2):", temp_E_R)
        E_R[i] = E_T0[i] - temp_E_R
        # println("radiated energy(MPa*m^2):",E_R[i])
        ER_M0_ratio[i] = -E_R[i]/(mu*area)    # without dimension along strike (in or out wall)
        # println("scaled energy:", ER_M0_ratio[i])

        seismic_moment[i] = mu*1e6*area*zdim    # assume that the rupture area is a square
        # println("seismic moment(N*m):", seismic_moment)
        # average stress drop
        del_sigma[i] = temp_sigma/zdim*3      # average stress drop for all asperities excluding the background matrix

        # radiation efficiency
        radiation_eff[i] = 2*mu/del_sigma[i]*ER_M0_ratio[i]      # here the unit of stress is MPa

        # average fault slip
        fault_slip[i] = area/zdim
        # rupture length along depth
        rupture_len[i] = -zdim

    end
    #  seismic_moment = filter!(x->x!=0, seismic_moment)
    #  del_sigma = filter!(x->x!=0, del_sigma)
    Mw = (2/3)*log10.(seismic_moment.*1e7) .- 10.7   # Kanamori(1977) and Hanks and Karamori(1979)

    return Mw, del_sigma, fault_slip, rupture_len, ER_M0_ratio, radiation_eff
end


function get_index(t, tStart, tEnd)     # get the index of time when earthquake begins and ends

    indx_start::Vector{Int64} = zeros(length(tStart[:,1]))
    indx_end::Vector{Int64} = zeros(length(tEnd[:,1]))
    println(typeof(t[1]))
    for i in eachindex(tStart[:,1])       # number of seismic events
        temp_start = findall(t[:] .<= tStart[i])[end]       # t[:] us recognized as float value but t is only a substring.
        temp_end = findall(t[:] .<= tEnd[i])[end]
        indx_start[i]= floor(temp_start/output_freq) + 1          # output every 10 timesteps
        indx_end[i] = floor(temp_end/output_freq)
    end
    return indx_start, indx_end

end

# seperate coseismic slip due to different events
function get_index_delfsec(N_events, delfsec)    # get the index of all coseismic slip when earthquake begins and ends
                                                # include the first artificial seismic event
                                  
    index_ds_start::Vector{Int} = zeros(N_events)      # 1d Vector
    index_ds_end::Vector{Int} = zeros(N_events)
    index_ds_start[1] = 1
    
    j = 1
    for i = 1:length(delfsec[:,1]) - 1
            if delfsec[i+1,1] - delfsec[i,1] >= 0.1    # the differnce of two group of 
                # coseismic slip is at least 1 m, so as to get the index of start and end of cosesimic slip
                index_ds_end[j] = i    
                index_ds_start[j+1] = i+1    
                j = j+1
            end
    end

    return index_ds_start, index_ds_end

end

# # get index of start of rupture() to plot the shear stress change
# function get_index(stress, taubefore, tauafter)        # timesteps, stress on points, 

#     len = length(taubefore[1,:])    # how many points on Fault
#     N = length(taubefore[:,1])    # how many seismic events on Fault
#     n = length(stress[:,1])      # how many timesteps/10
#     temp = 
#     for i in 1:N
#         temp_b = zeros(n)
#         temp_a = zeros(n)
#                 #  index_start[i] = findall(seismic_stress[len,:] .== taubefore[len,i])[1]
#         for j in 1:n        # calculate the differene between taubefore and shear stress at all timesteps/10
#             temp_b[j] = norm(stress[j,:] .- taubefore[i,:])
#             temp_a[j] = norm(stress[j,:] .- tauafter[i,:])
#         end

#         index_start[i] = findmin(temp_b)[2]
#         index_end[i] = findmin(temp_a)[2]
#     end

#     index_start, index_end
# end


# #.................................................
# # Compute the final Coseismic slip for each event
# #.................................................
# function Coslip(S, Slip, SlipVel, Stress, time_=zeros(1000000))
#     Vfmax = maximum(SlipVel, dims = 1)[:]

#     delfafter::Array{Float64,2} = zeros(size(Slip))
#     tStart::Array{Float64} = zeros(size(Slip[1,:]))
#     tEnd::Array{Float64} = zeros(size(Slip[1,:]))

#     taubefore::Array{Float64,2} = zeros(size(Slip))
#     tauafter::Array{Float64,2} = zeros(size(Slip))
    
#     hypo::Array{Float64} =  zeros(size(Slip[1,:]))   # Hypocenter
#     vhypo::Array{Float64} = zeros(size(Slip[1,:]))   # Velocity at hypocenter

#     Vthres = 0.001 # event threshold
#     slipstart = 0
#     it = 1; it2 = 1
#     delfref = zeros(size(Slip[:,1]))

#     for i = 1:length(Slip[1,:])

#         # Start of each event
#         if Vfmax[i] > 1.01*Vthres && slipstart == 0
#             delfref = Slip[:,i]
#             slipstart = 1
#             tStart[it2] = time_[i]
            
#             taubefore[:,it2] = Stress[:,i]
#             vhypo[it2], indx = findmax(SlipVel[:,i])

#             hypo[it2] = S.FltX[indx]

#             it2 = it2+1
#         end

#         # End of each event
#         if Vfmax[i] < 0.99*Vthres && slipstart == 1
#             delfafter[:,it] = Slip[:,i] - delfref
#             tauafter[:,it] = Stress[:,i]
#             tEnd[it] = time_[i]
#             slipstart = 0
#             it = it + 1
#         end
#     end

#     return delfafter[:,1:it-1], (taubefore-tauafter)[:,1:it-1], tStart[1:it2-1], tEnd[1:it-1], vhypo[1:it2-1], hypo[1:it2-1]
# end

