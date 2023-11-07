exp1(x::Float64) = ccall(:exp, Float64, (Float64,), x)
log1(x::Float64) = ccall(:log, Float64, (Float64,), x)

# IDstate functions
function IDS!(xLf, Vo, psi, dt, Vf, cnd, IDstate)   # default value is 2, aging law?

    #= compute slip-rates on fault based on different
       formulations =#
    if IDstate == 1
        psi1 = psi + dt*((Vo./xLf).*exp1(-psi) - abs(Vf)./xLf)     # lapusta (2000) for aging law: less accurate

        # psi1 = psi + dt* psi' = psi + dt*(theta'/theta) = psi + dt*(1/theta - Vf/xLf)

    # cnd = 1e-6
    # for V(t)dt/L < 1e-6 (low velocity), use the first three Taylor expansion of the second exponential 
    # for V(t)dt/L > 1e-6, use the full equation of theta(t+dt)

    elseif IDstate == 2      # aging law from kaneko(2011)
        VdtL = abs(Vf)*dt/xLf
        if VdtL < cnd        # linear form of the state variable evolution!
            psi1 = log1( exp1(psi-VdtL) + Vo*dt/xLf -
                        0.5*Vo*abs(Vf)*(dt^2)/(xLf^2))
        else
            psi1 = log1(exp1(psi-VdtL) + (Vo/abs(Vf))*(1-exp1(-VdtL)))
        end

    elseif IDstate == 3      # slip law from kaneko (2011)
        psi1 = exp1(-abs(Vf)*dt/xLf) * log1(abs(Vf)/Vo) + 
            exp1(-abs(Vf)*dt/xLf)*psi + log1(Vo/abs(Vf))

        # if ~any(imag(psi1)) == 0
        #     return abs(psi1)
        # end
    end

    return psi1

end

# On fault slip rates
# here we correct the slip rate from Vf1 to 0.5*(Vf1+Vf)
function IDS2!(xLf, Vo, psi, psi1, dt, Vf, Vf1, cnd, IDstate)
            
    if IDstate == 1
        psi2 = psi + 0.5*dt*( (Vo/xLf)*exp1(-psi) - abs(Vf)/xLf 
                                + (Vo/xLf)*exp1(-psi1) - abs(Vf1)/xLf )

    elseif IDstate == 2
        VdtL = 0.5*abs(Vf1 + Vf)*dt/xLf              # test the corrected sliprate: 0.5*abs(Vf1 + Vf)

        if VdtL < cnd
            psi2 = log1(exp1(psi-VdtL) + Vo*dt/xLf -
                            0.5*Vo*0.5*abs(Vf1 + Vf)*(dt)^2/(xLf^2))
  
        else
            psi2 = log1(exp1(psi-VdtL) + 
                            (Vo/(0.5*abs(Vf + Vf1)))*(1-exp1(-VdtL)))
        end

    elseif IDstate == 3
        psi2 = exp1(-0.5*abs(Vf + Vf1)*dt/xLf) * log1(0.5*abs(Vf + Vf1)/Vo) + 
                        exp1(-0.5*abs(Vf + Vf1)*dt/xLf)*psi + log1(Vo/(0.5*abs(Vf + Vf1)) )
    end

    return psi2
end

# Slip rates on fault for quasi-static regime
function slrFunc!(P::params_farray, NFBC, FltNglob, psi, psi1, Vf, Vf1, IDstate, tau1, dt, η)
    
    tauAB::Vector{Float64} = zeros(FltNglob)

    # for j = NFBC:FltNglob
    for j = NFBC[1]: NFBC[2]

        #  infer the next state variable from the previous sliprate Vf!
        psi1[j] = IDS!(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-6, IDstate)

        # tauAB[j] = tau1[j] + P.tauo[j]    # total shear stress 
        tauAB[j] = tau1[j] + P.tauo[j] - (η/2)*Vf[j]    # total shear stress - radiation damping      quasi-dynamic approximation

        # compute sliprate from a friction law together with total shear stress 
        fa = tauAB[j]/(P.Seff[j]*P.cca[j])      # tau_total/(a*sigma)

        help = -(P.fo[j] + P.ccb[j]*psi1[j])/P.cca[j]            # psi1 should not be changed!!

        # help = -(P.fo[j] + P.ccb[j]*psi1[j])/P.cca[j]
        help1 = exp(help + fa)
        help2 = exp(help - fa)
        # update the sliprate
        Vf1[j] = P.Vo[j]*(help1 - help2)
    end

return psi1, Vf1

end
