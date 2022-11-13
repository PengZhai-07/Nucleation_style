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
function moment_magnitude_new(mu, FltX, delfafter, stressdrops)
    # Final coseismic slip of each earthquake
    FltNglob = length(FltX)         # number of GLL nodes on fault

    iter = length(delfafter[1,:])    # number of events
    seismic_moment = zeros(iter)
    rupture_len = zeros(iter)
    fault_slip = zeros(iter)
    temp_sigma = 0
    iter2 = 1 

    del_sigma = zeros(iter)
    
    # calculate the nearest grid size
    dx = diff(FltX).*1e3        # Fltx: real depth(km) of all nodes on the fault 

    for i = 1:iter
        
        # slip threshold = 1% of maximum slip
        slip_thres = 0.01*maximum(delfafter[:,i])

        # area = slip*(rupture dimension along depth)
        # zdim = rupture along z dimension = depth rupture dimension
        area = 0; zdim = 0; temp_sigma = 0; temp_slip = 0

        for j = 1:FltNglob
            if delfafter[j,i] >= slip_thres
                area = area + delfafter[j,i]*dx[j-1]
                zdim = zdim + dx[j-1]     
                temp_slip = temp_slip + delfafter[j,i]*dx[j-1]

                # Avg. stress drops along rupture area
                temp_sigma = temp_sigma + stressdrops[j,i]*dx[j-1]
            end
        end
        
        seismic_moment[i] = mu*area*zdim
        # average stress drop
        del_sigma[i] = temp_sigma/zdim
        # average fault slip
        fault_slip[i] = temp_slip/zdim
        # rupture length along depth
        rupture_len[i] = -zdim

    end
    #  seismic_moment = filter!(x->x!=0, seismic_moment)
    #  del_sigma = filter!(x->x!=0, del_sigma)
    Mw = (2/3)*log10.(seismic_moment.*1e7) .- 10.7   # Kanamori(1977) and Hanks and Karamori(1979)

    return Mw, del_sigma, fault_slip, rupture_len
end


function get_index(t, tStart, tEnd)     # get the index of time when earthquake begins and ends

    indx_start::Vector{Int64} = zeros(length(tStart[:,1]))
    indx_end::Vector{Int64} = zeros(length(tEnd[:,1]))

    for i = 1 : length(tStart[:,1])       # number of seismic events
        temp_start = findall(t .<= tStart[i])[end]   
        temp_end = findall(t .<= tEnd[i])[end]
        indx_start[i]= floor(temp_start/10)
        indx_end[i] = floor(temp_end/10)
    end
    return indx_start, indx_end

end

# seperate coseismic slip due to different events
function get_index_delfsec(N_events, delfsec)    # get the index of coseismic slip when earthquake begins and ends

    index_ds_start::Vector{Int} = zeros(N_events)
    index_ds_end::Vector{Int} = zeros(N_events)
    index_ds_start[1] = 1
    index_ds_end[end] = size(delfsec)[1]
    j = 1
    for i = 1:length(delfsec[:,1])-1
            if delfsec[i+1,1] - delfsec[i,1] >= 0.1
                    index_ds_start[j+1] = i+1
                    index_ds_end[j] = i
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

