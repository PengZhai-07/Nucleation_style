####################################
#   fault slip behavior
####################################

# Fault Boundary function
function FBC!(IDstate::Int64, P::params_farray, NFBC, FltNglob, psi1, Vf1, tau1, psi2, Vf2, tau2, psi, Vf, FltVfree, dt)

    #  tauNR::Vector{BigFloat} = zeros(FltNglob)
    tauNR::BigFloat = 0.

    for j = NFBC[1]: NFBC[2]

        tauNR = 0.
        psi1[j] = IDS!(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-6, IDstate)          # state variable evolution

        Vf1[j], tau1[j] = NRsearch!(P.fo[j], P.Vo[j], P.cca[j], P.ccb[j], P.Seff[j],
                                    tauNR, P.tauo[j], psi1[j], P.FltZ[j], FltVfree[j])
    
        if Vf1[j] > 1e10 || isnan(Vf1[j]) == 1 || isnan(tau1[j]) == 1         
            
            println(" Fault Location = ", j)
            println(" Vf1 = ", Vf1[j])
            println(" tau1 = ", tau1[j])

            println("psi =", psi[j])
            println("psi1 =", psi1[j])
            # Save simulation results
            #filename = string(dir, "/data", name, "nrfail.jld2")
            #@save filename results(Stress,SlipVel, Slip, time_) 
            @error("NR SEARCH FAILED!")
            return
        end
        
        # correct the state variable again!!
        psi2[j] = IDS2!(P.xLf[j], P.Vo[j], psi[j], psi1[j], dt, Vf[j], Vf1[j], 1e-6, IDstate)
        
        # NRsearch 2nd loop: using new state variable psi2
        Vf2[j], tau2[j] = NRsearch!(P.fo[j], P.Vo[j], P.cca[j], P.ccb[j], P.Seff[j],
                                  tau1[j], P.tauo[j], psi2[j], P.FltZ[j], FltVfree[j])

    end

    return psi1, Vf1, tau1, psi2, Vf2, tau2
end

####################################
#   NEWTON RHAPSON SEARCH METHOD
####################################
# Newton Rhapson search method: find the sliprate and state shear stress satisfying a friction law 
function NRsearch!(fo, Vo, cca, ccb, Seff, tau, tauo, psi, FltZ, FltVfree)

    # psi maybe psi1 or psi2

    Vw = 1e10       
    
    fact = 1. + (Vo/Vw)*exp(-psi)      # fact = 1 + L/(Vw*theta)    # this is a number which is very close to 1
    fa::BigFloat = 0.
    help1::BigFloat = 0.
    help2::BigFloat = 0.
    delta::BigFloat = 0.

    # NR search point by point for tau if Vf < Vlimit

    # eps = 0.001*cca*Seff
    eps = 1e-6*cca*Seff
    k = 0
    delta = Inf

    while abs(delta) > eps

        fa = fact*tau/(Seff*cca)
        help = -(fo + ccb*psi)/cca

        help1 = exp(help + fa)
        help2 = exp(help - fa)

        Vf = Vo*(help1 - help2)         # calculate the slip rate first time

        # Vfprime is the first derivative of Vf to the traction(shear stress tau)
        Vfprime = fact*(Vo/(cca*Seff))*(help1 + help2)      

        # Our target loss function is (FltZ*FltVfree - FltZ*Vf + tauo - tau), 
        # which is a function of shear stress(traction): tau
        # residual shear stress (will decrease as iteration time increases))
        delta = (FltZ*FltVfree - FltZ*Vf + tauo - tau)/(1 + FltZ*Vfprime)

        # (FltZ * FltVfree + tauo - tau) should be zero, then 

        tau = tau + delta      # correct the shear stress
        k = k + 1       # record the number of searching

        if abs(delta) > 1e10 || k == 2000       # maximum search times
            println("k = ", k)
            # Save simulation results
            #filename = string(dir, "/data", name, "nrfail.jld2")
            #@save filename 
            #  @error("NR search fails to converge")

            return Float64(Vf), Float64(tau)
        end
    end
    
    # then we get the right shear stress and then we get the Vf
    fa = fact*tau/(Seff*cca)
    
    help = -(fo + ccb*psi)/cca

    help1 = exp(help + fa)
    help2 = exp(help - fa)

    Vf = Vo*(help1 - help2)

    return Float64(Vf), Float64(tau)
end
