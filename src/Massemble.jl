###################################################
#	ASSEMBLE THE MASS AND THE STIFFNESS MATRICES
###################################################

function Massemble!(NGLL, NelX, NelY, dxe, dye, ThickX, 
                   ThickY, rho1, vs1, rho2, vs2, iglob, 
                   M, x, y, jac)
    # Get the informaton about GLL points
    # May can be displaced by the input
    xgll, wgll, H = GetGLL(NGLL)
    wgll2 = wgll*wgll';

    rho::Matrix{Float64} = zeros(NGLL, NGLL)    # density 
    mu::Matrix{Float64} = zeros(NGLL, NGLL)     # shear modulus
    
    vso = zeros(NGLL, NGLL)
    vs = zeros(NGLL-1, NGLL)
    dx = zeros(NGLL-1, NGLL)
    muMax = 0
    dt = Inf
    
    #  # damage zone index
    #  damage_idx = zeros(Int, NelX*NelY)

    @inbounds @fastmath for ey = 1:NelY
        @inbounds @fastmath for ex = 1:NelX

            eo = (ey-1)*NelX + ex
            ig = iglob[:,:,eo]       # 2D index of No.eo element

            # Properties of heterogeneous medium
            if ex*dxe >= ThickX && (dye <= ey*dye <= ThickY)   #range of damage zone
                #  damage_idx[eo] = eo
                rho[:,:] .= rho2
                mu[:,:] .= rho2*vs2^2
            else          # range of host rock
                rho[:,:] .= rho1
                mu[:,:] .= rho1*vs1^2
            end
            # get the largest value of mu in all elements: the bigger one
            if  muMax < maximum(maximum(mu))
                muMax = maximum(maximum(mu))
            end

            # Diagonal Mass Matrix
            M[ig] .+= wgll2.*rho*jac

            # Local contributions to the stiffness matrix
            #  W[:,:,eo] .= wgll2.*mu;
            
            # Set timestep based on CFL criterion: dt<=dx/vs
            vso .= sqrt.(mu./rho)      # velocity of shear wave in all elements
            
            if dxe < dye   # row
                vs .= max.(vso[1:NGLL-1,:], vso[2:NGLL,:])    # largest row
                dx .= repeat( diff(xgll)*0.5*dxe, 1, NGLL)    # diff: get the difference : vso(i+1) - vso(i)
            else       # column
                vs .= max.(vso[:,1:NGLL-1], vso[:,2:NGLL])'    # largest column
                dx .= repeat( diff(xgll)*0.5*dye, 1, NGLL)
            end
            
            dtloc = dx./vs      # min(dx)/max(vs) : upper limit
            dt = minimum( push!(dtloc[1:end], dt) )    # minimum of all upper limits for all elements

        end
    end

    return M, dt, muMax 

end
