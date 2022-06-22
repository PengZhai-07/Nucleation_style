############################################
#   Compute the timestep for next iteration
############################################

function dtevol!(dt, dtmin, XiLf, FaultNglob, NFBC, Vf, isolver)
    # dt: current timstep
    dtmax::Int = 50 * 24 * 60*60		# 50 days
    dtincf::Float64 = 1.2

    if isolver == 1   # quasi-static regime

        # initial value of dt
        dtnx = dtmax

        # Adjust the timestep according to cell velocities and slip
        for i = NFBC:FaultNglob 

            if abs(Vf[i])*dtmax > XiLf[i]       # max slip in a timestep
                dtcell = XiLf[i]/abs(Vf[i])

                if dtcell < dtnx
                    dtnx = dtcell
                end
            end
        end

        if dtmin > dtnx
            dtnx = dtmin
        end
        # when to increase the timestep
        if dtnx > dtincf*dt
            dtnx = dtincf*dt
        end

        dt = dtnx

    elseif isolver == 2     # dynamic scheme
        
        dt = dtmin
    end

    return dt

end
