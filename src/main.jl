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
function healing2(t,tStart,dam, cos_reduction)
    """ hmax: coseismic damage amplitude
        r: healing rate (0.3454 => 20 years to heal completely)
                        (0.4605 => 15 years to heal completely)
                        (0.5756 => 12 years to heal completely)
                        (0.6908 => 10 years to heal completely)
                        (0.8635 => 8 years to healc completely)
                        (1.7269 => 4 years to heal completely)
                    """
    hmax = cos_reduction

    r =  0.8635   
    # t: current time of all simulation unit: seconds
    # tStart: time when earthquake happens  unit: seconds
    # dam : current ratio of damage zone and host rock(after coseismic rigidity reduction) 
    # 85% -> 80% -> 85%
    # when there is only 0.1% bias, we think the healing process finishes
    # ln(0.001)=-6.9078
    # healing time = 6.9078/r

    hmax*(1 .- exp.(-r*(t .- tStart)/P[1].yr2sec)) .+ dam     # hmax*(1 .- exp.(-r*(t .- tStart)/P[1].yr2sec)) > 0, when time is about 10 years, alphaa = dam + 0.05, healing completely!
end

function main(P, alphaa, cos_reduction, coseismic_b)
    # please refer to par.jl to see the specific meaning of P
    # P[1] = integer  Nel, FltNglob, yr2sec, Total_time, IDstate, nglob
    # P[2] = float    ETA, Vpl, Vthres, Vevne, dt
    # P[3] = float array   fo, Vo, xLf, M, BcBC, BcRC, FltL, FltZ, FltX, cca, ccb, Seff, tauo, XiLf, x_out, y_out
    # P[4] = integer array   iFlt, iBcB, iBcR, FltIglobBC, FltNI, out_seis
    # P[5] = ksparse   
    # P[6] = iglob, 
    # P[7] = NGLL
    # P[8] = wgll2
    # P[9] = nglob
    # P[10] = did

    #  ??? 
    #  W_orig = W[:,:,damage_idx]
    #  damage_amount::Float64 = 1.0

    # initial Shear modulus ratio of damage/host rock
    # pure elastic model: alphaa=1.0
    # immature fault: 80-85%
    # mature fault: 40%-45%

    # Time solver variables
    dt::Float64 = P[2].dt   # minimum timestep
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
    v .= 0.5e-3         #  half intial velocity on whole model: the real velocity is 1e-3 m/s
    # this initial velocity is only good for 50 Mpa and less, for higher normal stress, need larger initial velocity
    a::Vector{Float64} = zeros(P[1].nglob)   # relation between fault stress and acceleration?

    #.....................................
    # Stresses and time related variables on fault
    #.....................................
    
    FaultC::Vector{Float64} = zeros(P[1].FltNglob)     
    FltVfree::Vector{Float64} = zeros(length(P[4].iFlt))   # index of GLL nodes on the fault
    # velocity variables
    Vf::Vector{Float64} =  zeros(P[1].FltNglob)       
    Vf0::Vector{Float64} = zeros(length(P[4].iFlt))      # length(P[4].iFlt) = P[1].FltNglob
    Vf1::Vector{Float64} = zeros(P[1].FltNglob)
    Vf2::Vector{Float64} = zeros(P[1].FltNglob)

    # state variables
    psi::Vector{Float64} = zeros(P[1].FltNglob)
    #psi0::Vector{Float64} = zeros(P[1].FltNglob)       # intial state
    psi1::Vector{Float64} = zeros(P[1].FltNglob)
    psi2::Vector{Float64} = zeros(P[1].FltNglob)
    # stress variables
    tau::Vector{Float64} = zeros(P[1].FltNglob)
    tau1::Vector{Float64} = zeros(P[1].FltNglob)
    tau2::Vector{Float64} = zeros(P[1].FltNglob)
    tau3::Vector{Float64} = zeros(P[1].FltNglob)

    # Initial state variable
    # log.(P[3].Vo.*theta./xLf)

    # after earthquake, normal stress reduces and psi increase: how can normal stress influence state variable
    psi = P[3].tauo./(P[3].Seff.*P[3].ccb) - P[3].fo./P[3].ccb - (P[3].cca./P[3].ccb).*log.(2*v[P[4].iFlt]./P[3].Vo)
    #psi0 .= psi[:]
    # which kind of solver to use
    isolver::Int = 1         # quasi-static

    # # Some more initializations
    # r::Vector{Float64} = zeros(P[1].nglob)
    # beta_::Vector{Float64} = zeros(P[1].nglob)
    # alpha_::Vector{Float64} = zeros(P[1].nglob)

    # intial values on whole models!!!
    F::Vector{Float64} = zeros(P[1].nglob)
    dPre::Vector{Float64} = zeros(P[1].nglob)
    vPre::Vector{Float64} = zeros(P[1].nglob)
    dd::Vector{Float64} = zeros(P[1].nglob)
    dnew::Vector{Float64} = zeros(length(P[4].FltNI))  # save displacements for off-fault GLL nodes

    # Save output variables at certain timesteps: define output frequency
    tvsx::Float64 = 2e-0*P[1].yr2sec  # 2 years for interseismic period
    tvsxinc::Float64 = tvsx

    tevneinc::Float64 = 0.1    # 0.1 second for coseismic phase
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
    # Here v is not zero!!!  0.5e-3 m/s
    v = v[:] .- 0.5*P[2].Vpl   # initial slip rate on the whole model    ???
    Vf = 2*v[P[4].iFlt]      # about 1e-3
    iFBC::Vector{Int64} = findall(abs.(P[3].FltX) .> 24e3)   # index for points below the damage zone
    NFBC::Int64 = length(iFBC) + 1
    Vf[iFBC] .= 0.             # set the initial fault slip rate (within creeping fault) to be zero
    v[P[4].FltIglobBC] .= 0.   # set the initial fault slip rate (within creeping fault) to be zero
    # but intial slip rate on dynamic fault is not zero!!

    # on fault and off fault stiffness
    Ksparse = P[5]

    # Intact rock stiffness
    Korig = copy(Ksparse)   # K original： using with healing

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

    # Damage evolution stuff: using with healing 
    did = P[10]   # index of GLL nodes in fault damage zone
    dam = alphaa   # current rigidity ratio: initial value
    alpha_after = alphaa - cos_reduction
    print("alpha_after=", alpha_after,"\n")

    # Save parameters to file: from depth(48km) to shallow(0km)
    open(string(out_dir,"params.out"), "w") do io
        write(io, join(P[3].Seff/1e6, " "), "\n")  # unit: MPa
        write(io, join(P[3].tauo/1e6, " "), "\n")   # unit: MPa
        write(io, join(-P[3].FltX/1e3, " "), "\n")  # depth  unit: km   
        write(io, join(P[3].cca, " "), "\n")
        write(io, join(P[3].ccb, " "), "\n")
        write(io, join(P[3].xLf, " "), "\n")
    end
    open(string(out_dir,"mass_matrix.out"), "w") do io
        write(io, join(P[3].M, " "), "\n")  # unit: MPa
    end
    # Open files to begin writing
    open(string(out_dir,"stress.out"), "w") do stress    # shear stress 
    open(string(out_dir,"sliprate.out"), "w") do sliprate   # fault sliprate (Vpl+sliprate controlled by RSF)
    #open(string(out_dir,"slip.out"), "w") do slip   
    open(string(out_dir,"delfsec.out"), "w") do dfsec   # cultivate displacement(coseismic)
    open(string(out_dir,"delfsec_each_timestep.out"), "w") do dfsec_et   # cultivate displacement(coseismic)
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
    t = 0.  # current simualtion time
    Vfmax = 0.    # max slip rate of the fault
    
    # evolution of timesteps
    tStart2 = dt          # using with healing
    tStart = dt
    tEnd = dt
    # stress drop of earthquakes
    taubefore = P[3].tauo      
    tauafter = P[3].tauo
    # total dislocation(displacement) of dynamic fault after earthquake
    delfafter = 2*d[P[4].iFlt] .+ P[2].Vpl*t  
    hypo = 0.    # earthquake hypocenter
    d_hypo = 0.   # # cumulative slip at earthquake hypocenter
 
    # time dependence of b value
    SSS = 0
    seismogenic_depth = findall(abs(2e3) .< abs.(P[3].FltX) .<= abs(12e3))

    while t < P[1].Total_time
        it = it + 1

        t = t + dt   # dt is the initial smallest timestep: dtmin

        if isolver == 1    # quasi-static phase at the beginning!!!
            # intial velocity and displacement field for all GLL nodes
            vPre .= v     # non-zero
            dPre .= d     # Culmulative slip (zero) 

            Vf0 .= 2*v[P[4].iFlt] .+ P[2].Vpl   #  
            Vf  .= Vf0    # 0 m/s on kinematic fault segment while 1e-3 m/s on dynamic fault!

            # first two adjustation every time step during interseismic phase
            for p1 = 1:2
                
                # step 1: predict displacement on the fault
                # v[P[4].iFlt]: initial slip rate on fault, large velocity on kinematic fault 
                F .= 0.
                F[P[4].iFlt] .= dPre[P[4].iFlt] .+ v[P[4].iFlt]*dt

                # Assign previous displacement field as initial guess for off-fault nodes
                dnew .= d[P[4].FltNI]    # zero
                
                # step 2: solve for displacement in the medium
                # Solve d = K^-1 F by MGCG
                # stiffness matrix on the fault 
                rhs = (mul!(tmp,Ksparse,F))[P[4].FltNI]   # However comparing this to LinearAlgebra.mul! shows that the later is much faster.
                #  rhs = (Ksparse*F)[P[4].FltNI]

                # direct inversion
                #  dnew = -(kni\rhs)

                # mgcg
                dnew = cg!(dnew, kni, rhs, Pl=p, abstol=1e-6)    # note: in new version of julia, there is only abstol  
  
                # update displacement on the medium
                d[P[4].FltNI] .= dnew

                # make d = F on the fault: update displacement on fault
                d[P[4].iFlt] .= F[P[4].iFlt]

                # step3: Compute on-fault stress: f
                a .= 0.
                mul!(a,Ksparse,d)
                #   a = Ksparse*d

                # step4: compute traction on the fault
                # Enforce K*d to be zero for velocity boundary (0-20 km)
                # there is no traction on creeping fault 
                a[P[4].FltIglobBC] .= 0.       # for creeping fault, traction is zero 

                # step4
                tau1 .= -a[P[4].iFlt]./P[3].FltL
                
                # step5
                # Function to calculate on-fault sliprate on whole fault line

                # note that slrFunc! depends on ccb
                psi1, Vf1 = slrFunc!(P[3], NFBC, P[1].FltNglob, psi, psi1, Vf, Vf1, P[1].IDstate, tau1, dt)   # from other functions
                
                # step6: correct slip rate on the fault
                Vf1[iFBC] .= P[2].Vpl     # set slip rate on creep fault to be plate motion rate
                # current slip rate 
                Vf .= (Vf0 + Vf1)/2   # slip rate on creeping fault:  P[2].Vpl  
                v[P[4].iFlt] .= 0.5*(Vf .- P[2].Vpl)   # slip rate on creeping fault:  0 

            end
            
            # update current variable
            psi .= psi1[:]
            tau .= tau1[:]

            # # creeping fault
            tau[iFBC] .= 0.
            Vf[iFBC] .= P[2].Vpl

            # on fault GLL nodes
            v[P[4].iFlt] .= 0.5*(Vf .- P[2].Vpl)   

            # off-fault GLL nodes
            v[P[4].FltNI] .= (d[P[4].FltNI] .- dPre[P[4].FltNI])/dt

            # Line 731: P_MA: Omitted
            # reset the fault stress(or acceleration) to be zero
            a .= 0.
            # reset the culmulative displacement and velocity within creeping fault to be zero
            d[P[4].FltIglobBC] .= 0.
            v[P[4].FltIglobBC] .= 0.

            #---------------
            # Healing stuff: only for interseismic quasi-static phase
            # --------------
            # 

# comment this part for no healing: start  

            # if  it > 3
                
            #     #if t > 10*P[1].yr2sec     # healing after 10 year, neglect the first event
            #         alphaa = healing2(t, tStart2, dam, cos_reduction)           # healing from dam
            #         #  alphaa[it] = αD(t, tStart2, dam)
            #     #end

            #     for id in did
            #         Ksparse[id] = alphaa*Korig[id]   # define the stiffness of fault damage zone
            #     end

            #     #println("alpha healing = ", alphaa[it])

            #     # Linear solver stuff
            #     kni = -Ksparse[P[4].FltNI, P[4].FltNI]
            #     nKsparse = -Ksparse
            #     # multigrid
            #     ml = ruge_stuben(kni)
            #     p = aspreconditioner(ml)

            #     # faster matrix multiplication
            #     #  Ksparse = Ksparse'
            #     #  nKsparse = nKsparse'
            #     #  kni = kni'
            # end

# comment this part for no healing: end
        
        # If isolver != 1, or max slip rate is > 10^-3 m/s , dynamic phase
        else

            dPre .= d
            vPre .= v

            # step1: Update displacement and partial velocity
            d .= d .+ dt.*v .+ (half_dt_sq).*a

            # Prediction
            v .= v .+ half_dt.*a
            a .= 0.   # traction on fault

            # step2: computing the internal forces -K*d[t+1] stored in global array 'a'
            mul!(a,nKsparse,d)
            #   a = nKsparse*d

            # step3 computing the 'stick' traction
            # Enforce K*d to be zero for velocity boundary
            a[P[4].FltIglobBC] .= 0.     # creeping fault

            # Absorbing boundaries(Bottom and right)    ??? no stress or acceleration
            a[P[4].iBcB] .= a[P[4].iBcB] .- P[3].BcBC.*v[P[4].iBcB]
            a[P[4].iBcR] .= a[P[4].iBcR] .- P[3].BcRC.*v[P[4].iBcR]

            ###### Fault Boundary Condition: Rate and State #############
            FltVfree .= 2*v[P[4].iFlt] .+ 2*half_dt*a[P[4].iFlt]./P[3].M[P[4].iFlt]
            Vf .= 2*vPre[P[4].iFlt] .+ P[2].Vpl

            # step4: find fault traction and slip velocity satisfying a friction law and the relation
            # Sliprate and NR search
            psi1, Vf1, tau1, psi2, Vf2, tau2 = FBC!(P[1].IDstate, P[3], NFBC, P[1].FltNglob, psi1, Vf1, tau1, psi2, Vf2, tau2, psi, Vf, FltVfree, dt)

            # step5: add the fault boundary term to the sum of internal forces
            tau .= tau2 .- P[3].tauo
            tau[iFBC] .= 0.
            psi .= psi2
            # 
            a[P[4].iFlt] .= a[P[4].iFlt] .- P[3].FltL.*tau
            ########## End of fault boundary condition ##############


            # step6: Solve for a_new acceleration 
            a .= a./P[3].M

            # step7: Correction
            v .= v .+ half_dt*a
            v .= v .- 0.5*P[2].Vpl       # get the velocity perturbation

            v[P[4].FltIglobBC] .= 0.
            a[P[4].FltIglobBC] .= 0.

            #### Line 861: Omitting P_Ma


        end # of isolver if loop

        Vfmax = 2*maximum(v[P[4].iFlt]) .+ P[2].Vpl   # background plate motion rate: P[2].Vpl

# velocity dependent b (evolution effect)

        if t > 10*P[1].yr2sec

            if  Vfmax <= 1e-5          
                P[3].ccb[seismogenic_depth] .= 0.019

            elseif 1e-5 < Vfmax < 1e-3  
                # linear increase in Cartesian velocity
                # K_b = (0.025-0.019)/(1e-3-1e-5)
                # P[3].ccb[seismogenic_depth] .= 0.019 .+ (Vfmax - 1e-5)* K_b

                # sin increase in Cartesian velocity
                # normalization
                K_b = (1+sin((Vfmax - 1e-5)/(1e-3 - 1e-5)*pi-pi/2))/2
                P[3].ccb[seismogenic_depth] .= 0.019 + (coseismic_b -0.019) * K_b

                # linear increase in log scale velocity
                # K_b = (0.025-0.019)/(log10(1e-3) - log10(1e-5))
                # P[3].ccb[seismogenic_depth] .= 0.019 .+ (log10(Vfmax) - log(1e-5))* K_b  

                # if  SSS == 0
                #     println(P[3].ccb)
                # end
                # SSS = SSS + 1

            elseif Vfmax >= 1e-3
                P[3].ccb[seismogenic_depth] .= coseismic_b
            end    
        end

        #-----
        # Output the variables before and after events
        #-----
        # record the coseismic event!!
        if  Vfmax > 1.01*P[2].Vthres && slipstart == 0

            it_s = it_s + 1
            delfref = 2*d[P[4].iFlt] .+ P[2].Vpl*t
            
            slipstart = 1   # sign for new earthquake!!

            tStart = t
            taubefore = (tau +P[3].tauo)./1e6     # once switch into dynamic solution, record the shear stress before earthquake
            # hypocenter: fault location where slip rate exceed threshold value firstly!!
            vhypo, indx = findmax(2*v[P[4].iFlt] .+ P[2].Vpl)
            hypo = P[3].FltX[indx]
            d_hypo = delfref[indx]

        end

        # slip rate is lower than Vthres(When earthquake ends or at the beginning )
        if Vfmax < 0.99*P[2].Vthres && slipstart == 1
            it_e = it_e + 1
            delfafter = 2*d[P[4].iFlt] .+ P[2].Vpl*t .- delfref     # coseismic slip
            
            tEnd = t 
            tauafter = (tau +P[3].tauo)./1e6    # once switch into quasi-static solution, record the shear stress
            
            # Save start and end time and stress
            write(event_time, join(hcat(tStart,tEnd, -hypo, d_hypo), " "), "\n")
            write(event_stress, join(hcat(taubefore, tauafter), " "), "\n")
            write(dfafter, join(delfafter, " "), "\n")
            
            slipstart = 0
            
# note comment this part without healing: start   
                # # at the end of each earthquake, the shear wave velocity in the damaged zone reduces

                # # Time condition of 10 years
                # #if t > 10*P[1].yr2sec 
                #     #  use this for no permanent damage    65-60-65%
                #        alphaa = alpha_after   # a constant value after event
                #        dam = alphaa    # current rigidity ratio use for healing

                #     #  Use this for permanent damage: 1%
                #     #  alphaa = alphaa - 0.06      
                #     #  dam = alphaa
                #     #  if dam < 0.60     # lowest rigidity ratio
                #         #  alphaa = 0.60
                #         #  dam = 0.60
                #     #  end
                
                # # it's necessary to change the stiffness matrix if the rigidity ratio is changed 
                #     tStart2 = t            # used for healing!

                #     for id in did
                #         Ksparse[id] = alphaa*Korig[id]      # calculate the stiffness of fault damage zone again
                #     end

                #     # Linear solver stuff
                #     kni = -Ksparse[P[4].FltNI, P[4].FltNI]
                #     nKsparse = -Ksparse
                #     # multigrid
                #     ml = ruge_stuben(kni)
                #     p = aspreconditioner(ml)

                # #end

# note comment this part without healing: end  

                println("alphaa = ", alphaa)   # output the rigidity ratio after every earthquake 
                # warning: println is different from @printf !!!
        end
        
        #-----
        # Output the variables certain timesteps: 2yr interseismic, 1 sec coseismic
        #-----
        if t > tvsx                  # 2years
            ntvsx = ntvsx + 1
            idd += 1
            #  write(stress, join((tau + P[3].tauo)./1e6, " "), "\n")
            write(dfyr, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")

            tvsx = tvsx + tvsxinc     # output frequency: 2 years
        end

        if Vfmax > 1.01*P[2].Vthres

            if mod(it,10) == 0             # output the shear stress every 10 steps as the shear stess
                write(dfsec_et, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")
            end

            if idelevne == 0                  # record the first step
                nevne = nevne + 1
                idd += 1
                idelevne = 1
                tevneb = t
                tevne = tevneinc      # 0.1s

                write(dfsec, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")
            
            end

            if idelevne == 1 && (t - tevneb) > tevne      # output at least every 0.1 s in the following steps
                nevne = nevne + 1
                idd += 1

                write(dfsec, join(2*d[P[4].iFlt] .+ P[2].Vpl*t, " "), "\n")
                tevne = tevne + tevneinc
            end

        else
            idelevne = 0
        end

        current_sliprate = 2*v[P[4].iFlt] .+ P[2].Vpl

        # Output timestep info on screen: every 50 timesteps
        if mod(it,500) == 0
            @printf("Time (yr) = %1.5g\n", t/P[1].yr2sec) 
            #  println("Vfmax = ", maximum(current_sliprate))
        end

        # Write stress, sliprate, slip to file every 5 timesteps
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
        write(Vf_time, join(hcat(t, Vfmax, Vf[end], alphaa, isolver, maximum(P[3].ccb)), " "), "\n")
        
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
end

# These ends are for opened output files