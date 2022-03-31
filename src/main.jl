###############################################################################
#
#	SPECTRAL ELEMENT METHOD FOR EARTHQUAKE CYCLE SIMULATION
#
#   Written in: Julia 1.0
#
#	Created: 09/18/2020
#   Author: Prithvi Thakur (Original code by Kaneko et al. 2011)
#
#	Adapted from Kaneko et al. (2011)
#	and J.P. Ampuero's SEMLAB
#
###############################################################################

# Healing exponential function
function healing2(t,tStart,dam)
    """ hmax: coseismic damage amplitude
        r: healing rate (0.05 => 138 years to heal completely)
                        (0.5 => 13.81 years to heal completely)
                        (0.7 => 9.87 years to heal completely)
                        (0.8 => 8.63 years to heal completely)
    """
    hmax = 0.05
    r =  0.7   # 1/1.5
    # t: current time of all simulation unit: seconds
    # tStart: time when earthquake happens  unit: seconds
    # dam : current ratio of damage zone and host rock(after coseismic rigidity reduction) 
    # 65% -> 60% -> 65%
    # when there is only 0.1% bias, we think the healing process finishes
    # log(0.001)=-6.9078
    hmax*(1 .- exp.(-r*(t .- tStart)/P[1].yr2sec)) .+ dam     # hmax*(1 .- exp.(-r*(t .- tStart)/P[1].yr2sec)) > 0, when time is 10years, alphaa = dam + 0.05
end

function main(P,alphaa)
    # please refer to par.jl to see the specific meaning of P
    # P[1] = integer  Nel, FltNglob, yr2sec, Total_time, IDstate, nglob
    # P[2] = float    ETA, Vpl, Vthres, Vevne, dt
    # P[3] = float array   fo, Vo, xLf, M, BcBC, BcRC, FltL, FltZ, FltX, cca, ccb, Seff, tauo, XiLf, x_out, y_out
    # P[4] = integer array   iFlt, iBcB, iBcR, FltIglobBC, FltNI, out_seis
    # P[5] = ksparse   

    #  ??? 
    #  W_orig = W[:,:,damage_idx]
    #  damage_amount::Float64 = 1.0

    # initial Shear modulus ratio of damage/host rock
    # pure elastic model: alphaa=1.0
    # immature fault: 80-85%
    # mature fault: 40%-45%

    # Time solver variables
    dt::Float64 = P[2].dt   # dt is variable at different time of earthquake cycle
    dtmin::Float64 = dt    
    half_dt::Float64 = 0.5*dtmin
    half_dt_sq::Float64 = 0.5*dtmin^2

    # dt modified slightly for damping (if with viscosity, the timestep is smaller)
    if P[2].ETA != 0
        dt = dt/sqrt(1 + 2*P[2].ETA)
    end

    # Initialize kinematic field: global arrays
    d::Vector{Float64} = zeros(P[1].nglob)   # initial displacement
    v::Vector{Float64} = zeros(P[1].nglob)
    v .= 0.5e-3         # half of Vthres(1e-3 m/second)   intial velocity(but not enough to lead to a earthquake)
    a::Vector{Float64} = zeros(P[1].nglob)  #???

    #.....................................
    # Stresses and time related variables on fault
    #.....................................
    # specific meanings
    
    FaultC::Vector{Float64} = zeros(P[1].FltNglob)
    FltVfree::Vector{Float64} = zeros(length(P[4].iFlt))
    Vf::Vector{Float64} =  zeros(P[1].FltNglob)
    Vf0::Vector{Float64} = zeros(length(P[4].iFlt))   # length(P[4].iFlt) = P[1].FltNglob
    Vf1::Vector{Float64} = zeros(P[1].FltNglob)
    Vf2::Vector{Float64} = zeros(P[1].FltNglob)
    # state variables
    psi::Vector{Float64} = zeros(P[1].FltNglob)
    psi0::Vector{Float64} = zeros(P[1].FltNglob)   # intial state
    psi1::Vector{Float64} = zeros(P[1].FltNglob)
    psi2::Vector{Float64} = zeros(P[1].FltNglob)
    # stress variables
    tau::Vector{Float64} = zeros(P[1].FltNglob)
    tau1::Vector{Float64} = zeros(P[1].FltNglob)
    tau2::Vector{Float64} = zeros(P[1].FltNglob)
    tau3::Vector{Float64} = zeros(P[1].FltNglob)

    # Initial state variable
    psi = P[3].tauo./(P[3].Seff.*P[3].ccb) - P[3].fo./P[3].ccb - (P[3].cca./P[3].ccb).*log.(2*v[P[4].iFlt]./P[3].Vo)
    psi0 .= psi[:]
  # which kind of solver to use
    isolver::Int = 1

    # # Some more initializations
    # r::Vector{Float64} = zeros(P[1].nglob)
    # beta_::Vector{Float64} = zeros(P[1].nglob)
    # alpha_::Vector{Float64} = zeros(P[1].nglob)

    F::Vector{Float64} = zeros(P[1].nglob)
    dPre::Vector{Float64} = zeros(P[1].nglob)
    vPre::Vector{Float64} = zeros(P[1].nglob)
    dd::Vector{Float64} = zeros(P[1].nglob)
    dnew::Vector{Float64} = zeros(length(P[4].FltNI))  # off-fault GLL nodes

    # Save output variables at certain timesteps: define output frequency
    tvsx::Float64 = 2e-0*P[1].yr2sec  # 2 years for interseismic period
    tvsxinc::Float64 = tvsx

    tevneinc::Float64 = 0.1    # 0.5 second for seismic period
    delfref = zeros(P[1].FltNglob)

    # Iterators
    idelevne::Int= 0
    tevneb::Float64= 0.
    tevne::Float64= 0.
    ntvsx::Int= 0
    nevne::Int= 0
    slipstart::Int= 0
    idd::Int = 0
    it_s = 0; it_e = 0
    rit = 0

    v = v[:] .- 0.5*P[2].Vpl   # initial slip rate
    # Vf = 2*v[P[4].iFlt]  # Plate motion rate?
    iFBC::Vector{Int64} = findall(abs.(P[3].FltX) .> 24e3)   # index for points below the damage zone
    NFBC::Int64 = length(iFBC) + 1
    Vf[iFBC] .= 0.   # set the initial fault slip rate(below the damage zone) to be zero

    v[P[4].FltIglobBC] .= 0.   # Reset the initial velocity to be zero ???

    # on fault and off fault stiffness
    Ksparse = P[5]

    # Intact rock stiffness
    Korig = copy(Ksparse)   # K original

    # Linear solver stuff
    kni = -Ksparse[P[4].FltNI, P[4].FltNI]   # stiffness of off-fault GLL nodes
    nKsparse = -Ksparse
    
    # algebraic multigrid preconditioner
    ml = ruge_stuben(kni)
    p = aspreconditioner(ml)
    tmp = copy(a)

    # faster matrix multiplication
     #  Ksparse = Ksparse'
     #  nKsparse = nKsparse'
     #  kni = kni'

    #  Ksparse = ThreadedMul(Ksparse)
    #  nKsparse = ThreadedMul(nKsparse)
    #  kni = ThreadedMul(kni)

    # Damage evolution stuff
    did = P[10]   # index of GLL nodes in fault damage zone
    dam = alphaa   # initial damage ratio

    # Save parameters to file: from depth(48km) to shallow(0km)
    open(string(out_dir,"params.out"), "w") do io
        write(io, join(P[3].Seff/1e6, " "), "\n")  # unit: MPa
        write(io, join(P[3].tauo/1e6, " "), "\n")   # unit: MPa
        write(io, join(-P[3].FltX/1e3, " "), "\n")  # depth  unit: km   
        write(io, join(P[3].cca, " "), "\n")
        write(io, join(P[3].ccb, " "), "\n")
        write(io, join(P[3].xLf, " "), "\n")
    end

    # Open files to begin writing
    open(string(out_dir,"stress.out"), "w") do stress    # shear stress 
    open(string(out_dir,"sliprate.out"), "w") do sliprate   # fault sliprate (Vpl+sliprate controlled by RSF)
    #open(string(out_dir,"slip.out"), "w") do slip   
    open(string(out_dir,"delfsec.out"), "w") do dfsec   # cultivate displacement(coseismic)
    open(string(out_dir,"delfyr.out"), "w") do dfyr   # cultivate displacement(interseismic)
    open(string(out_dir,"event_time.out"), "w") do event_time    # start and end time and hypocenter of earthquake event
    open(string(out_dir,"event_stress.out"), "w") do event_stress  # shear stress before and after 
    open(string(out_dir,"coseismic_slip.out"), "w") do dfafter
        # time/max slip rate on fault/slip rate near ground/rigidity ratio
    open(string(out_dir,"time_velocity.out"), "w") do Vf_time  

    #....................
    # Start of time loop
    #....................
    it = 0  # current the number of timesteps
    t = 0.  # current time
    Vfmax = 0.    # max slip rate of the fault
    
    tStart2 = dt
    tStart = dt
    tEnd = dt
    # stress drop before and after earthquake
    taubefore = P[3].tauo      
    tauafter = P[3].tauo
    # background loading rate
    delfafter = 2*d[P[4].iFlt] .+ P[2].Vpl*t  
    hypo = 0.    # earthquake location

    while t < P[1].Total_time
        it = it + 1
        t = t + dt   # dt is the initial smallest timestep
        
        if isolver == 1    # quasi-static phase at the beginning!!!
            vPre .= v
            dPre .= d

            Vf0 .= 2*v[P[4].iFlt] .+ P[2].Vpl   # initial fault slip rate
            Vf  .= Vf0    # 1e-3

            # first two adjustation every time step during interseismic phase
            for p1 = 1:2

                # Compute the on-Fault displacement
                F .= 0.
                F[P[4].iFlt] .= dPre[P[4].iFlt] .+ v[P[4].iFlt]*dt

                # Assign previous displacement field as initial guess
                dnew .= d[P[4].FltNI]


                # Solve d = K^-1F by MGCG
                rhs = (mul!(tmp,Ksparse,F))[P[4].FltNI]
                #  rhs = (Ksparse*F)[P[4].FltNI]

                # direct inversion
                #  dnew = -(kni\rhs)

                # mgcg
                dnew = cg!(dnew, kni, rhs, Pl=p, abstol=1e-6)    # note: in new version of julia, there is only abstol  

                # update displacement on the medium
                d[P[4].FltNI] .= dnew

                # make d = F on the fault
                d[P[4].iFlt] .= F[P[4].iFlt]

                # Compute on-fault stress
                a .= 0.
                mul!(a,Ksparse,d)
                #   a = Ksparse*d

                # Enforce K*d to be zero for velocity boundary
                a[P[4].FltIglobBC] .= 0.

                tau1 .= -a[P[4].iFlt]./P[3].FltL

                # Function to calculate on-fault sliprate
                psi1, Vf1 = slrFunc!(P[3], NFBC, P[1].FltNglob, psi, psi1, Vf, Vf1, P[1].IDstate, tau1, dt)   # from other functions
                
                # set the initial fault slip rate(below the damage zone) to be plate motion rate
                Vf1[iFBC] .= P[2].Vpl     # slip rate on kinematic fault
                Vf .= (Vf0 + Vf1)/2   # >24km: v[P[4].iFlt] .+ P[2].Vpl  
                v[P[4].iFlt] .= 0.5*(Vf .- P[2].Vpl)  # 0.5 * v[P[4].iFlt]

            end

            psi .= psi1[:]
            tau .= tau1[:]
            # kinematic fault
            tau[iFBC] .= 0.
            Vf1[iFBC] .= P[2].Vpl
            # on fault GLL nodes
            v[P[4].iFlt] .= 0.5*(Vf1 .- P[2].Vpl)   # keep the slip rate of kinematic fault as a constant!!
            # off-fault GLL nodes
            v[P[4].FltNI] .= (d[P[4].FltNI] .- dPre[P[4].FltNI])/dt

            # Line 731: P_MA: Omitted
            # specific meaning??
            a .= 0.
            d[P[4].FltIglobBC] .= 0.
            v[P[4].FltIglobBC] .= 0.

            #---------------
            # Healing stuff: Ignore for now
            # --------------
            # Normal calculation from 3rd step!!!
            if  it > 3
                #if t > 10*P[1].yr2sec     # healing after 10 year, neglect the first event
                    #alphaa = healing2(t, tStart2, dam)
                    #  alphaa[it] = αD(t, tStart2, dam)
                #end

                for id in did
                    Ksparse[id] = alphaa*Korig[id]   # define the stiffness of fault damage zone
                end

                #println("alpha healing = ", alphaa[it])

                # Linear solver stuff
                kni = -Ksparse[P[4].FltNI, P[4].FltNI]
                nKsparse = -Ksparse
                # multigrid
                ml = ruge_stuben(kni)
                p = aspreconditioner(ml)

                # faster matrix multiplication
                #  Ksparse = Ksparse'
                #  nKsparse = nKsparse'
                #  kni = kni'
            end

        
        # If isolver != 1, or max slip rate is > 10^-3 m/s , dynamic phase
        else

            dPre .= d
            vPre .= v

            # Update
            d .= d .+ dt.*v .+ (half_dt_sq).*a

            # Prediction
            v .= v .+ half_dt.*a
            a .= 0.

            # Internal forces -K*d[t+1] stored in global array 'a'
            mul!(a,nKsparse,d)
            #   a = nKsparse*d

            # Enforce K*d to be zero for velocity boundary
            a[P[4].FltIglobBC] .= 0.

            # Absorbing boundaries(Bottom and right)
            a[P[4].iBcB] .= a[P[4].iBcB] .- P[3].BcBC.*v[P[4].iBcB]
            a[P[4].iBcR] .= a[P[4].iBcR] .- P[3].BcRC.*v[P[4].iBcR]

            ###### Fault Boundary Condition: Rate and State #############
            FltVfree .= 2*v[P[4].iFlt] .+ 2*half_dt*a[P[4].iFlt]./P[3].M[P[4].iFlt]
            Vf .= 2*vPre[P[4].iFlt] .+ P[2].Vpl


            # Sliprate and NR search
            psi1, Vf1, tau1, psi2, Vf2, tau2 = FBC!(P[1].IDstate, P[3], NFBC, P[1].FltNglob, psi1, Vf1, tau1, psi2, Vf2, tau2, psi, Vf, FltVfree, dt)

            tau .= tau2 .- P[3].tauo
            tau[iFBC] .= 0.
            psi .= psi2
            a[P[4].iFlt] .= a[P[4].iFlt] .- P[3].FltL.*tau
            ########## End of fault boundary condition ##############


            # Solve for a_new
            a .= a./P[3].M

            # Correction
            v .= v .+ half_dt*a

            v[P[4].FltIglobBC] .= 0.
            a[P[4].FltIglobBC] .= 0.

            #### Line 861: Omitting P_Ma


        end # of isolver if loop

        Vfmax = 2*maximum(v[P[4].iFlt]) .+ P[2].Vpl   # background plate motion rate: P[2].Vpl

        #-----
        # Output the variables before and after events
        #-----
        # record the coseismic event!!
        if  Vfmax > 1.01*P[2].Vthres && slipstart == 0
            it_s = it_s + 1
            delfref = 2*d[P[4].iFlt] .+ P[2].Vpl*t
            
            slipstart = 1   # sign for new earthquake!!

            tStart = t
            taubefore = (tau +P[3].tauo)./1e6
            # hypocenter: fault location where slip rate exceed threshold value firstly!!
            vhypo, indx = findmax(2*v[P[4].iFlt] .+ P[2].Vpl)
            hypo = P[3].FltX[indx]

        end

        # slip rate is lower than Vthres(When earthquake ends or at the beginning )
        if Vfmax < 0.99*P[2].Vthres && slipstart == 1
            it_e = it_e + 1
            delfafter = 2*d[P[4].iFlt] .+ P[2].Vpl*t .- delfref     # coseismic slip
            
            tEnd = t 
            tauafter = (tau +P[3].tauo)./1e6
            
            # Save start and end time and stress
            write(event_time, join(hcat(tStart,tEnd, -hypo), " "), "\n")
            write(event_stress, join(hcat(taubefore, tauafter), " "), "\n")
            write(dfafter, join(delfafter, " "), "\n")
            
            slipstart = 0
            
            # # at the end of each earthquake, the shear wave velocity in the damaged zone reduces by 5%

            #     # Time condition of 10 years
            #     #if t > 10*P[1].yr2sec 
            #         #  use this for no permanent damage    65-60-65%
            #             alphaa = 0.60
            #             dam = alphaa      # rigidity ratio before healing 


            #         #  Use this for permanent damage
            #         #  alphaa = alphaa - 0.05
            #         #  dam = alphaa
            #         #  if dam < 0.60     # lowest rigidity ratio
            #             #  alphaa = 0.60
            #             #  dam = 0.60
            #         #  end
                
            #     # it's necessary to change the stiffness matrix if the rigidity ratio is changed 
            #         tStart2 = t            # used for healing!

            #         for id in did
            #             Ksparse[id] = alphaa*Korig[id]      # calculate the stiffness of fault damage zone again
            #         end

            #         # Linear solver stuff
            #         kni = -Ksparse[P[4].FltNI, P[4].FltNI]
            #         nKsparse = -Ksparse
            #         # multigrid
            #         ml = ruge_stuben(kni)
            #         p = aspreconditioner(ml)

                #end

                println("alphaa = ", alphaa)   # output the rigidity ratio after every earthquake 

            #  end

        end
        
        #-----
        # Output the variables certain timesteps: 2yr interseismic, 1 sec coseismic
        #-----
        if t > tvsx
            ntvsx = ntvsx + 1
            idd += 1
            #  write(stress, join((tau + P[3].tauo)./1e6, " "), "\n")
            write(dfyr, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")

            tvsx = tvsx + tvsxinc
        end

        if Vfmax > P[2].Vevne
            if idelevne == 0
                nevne = nevne + 1
                idd += 1
                idelevne = 1
                tevneb = t
                tevne = tevneinc

                #  write(stress, join((tau + P[3].tauo)./1e6, " "), "\n")
                write(dfsec, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")
            end

            if idelevne == 1 && (t - tevneb) > tevne
                nevne = nevne + 1
                idd += 1

                write(dfsec, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")
                tevne = tevne + tevneinc
            end

        else
            idelevne = 0
        end

        current_sliprate = 2*v[P[4].iFlt] .+ P[2].Vpl

        # Output timestep info on screen: every 500 timesteps
        if mod(it,500) == 0
            @printf("Time (yr) = %1.5g\n", t/P[1].yr2sec) 
            #  println("Vfmax = ", maximum(current_sliprate))
        end


        # Write stress, sliprate, slip to file every 10 timesteps
        if mod(it,10) == 0
            write(sliprate, join(2*v[P[4].iFlt] .+ P[2].Vpl, " "), "\n")
            write(stress, join((tau + P[3].tauo)./1e6, " "), "\n")
        end

        # Determine quasi-static or dynamic regime based on max-slip velocity
        #  if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3
        # when to change the solver
        if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3    
            # 0.5e-3 is the initial slip rate, so that there is an initial earthquake at zero time!!
            # in addition, 5e-3 is half of the vthres, if it necessary to convert to dynamic regime in advance??
            isolver = 1   # quasi-static
        else
            isolver = 2   # dynamic
        end

        # Write max sliprate and time
            # t: simulation time (seconds)
            # Vfmax: max slip rate on the fault
            # Vf[end] is fault slip rate on the surface!!
            # alphaa: current rigidity ratio of fault damage zone
        write(Vf_time, join(hcat(t, Vfmax, Vf[end], alphaa), " "), "\n")

        # Compute next timestep dt: adaptive!!
        dt = dtevol!(dt , dtmin, P[3].XiLf, P[1].FltNglob, NFBC, current_sliprate, isolver)


    end # end of time loop

    # close files
# end of writing data into files
end
end
end
end
end
end
end
end

end