# Material properties for a narrow rectangular damaged zone of half-thickness ThickY and depth ThickX 
# define the properties of fault damage zone using locations of elements: on the level of elements, no change of properties inside elements  
function MaterialProperties(NelX, NelY,NGLL, dxe, dye, ThickX, ThickY, wgll2, rho1, vs1, rho2, vs2, Domain)
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    W::Array{Float64,3} = zeros(NGLL, NGLL, NelX*NelY)
    
    # damage zone index
    #  damage_idx = zeros(Int, NelX*NelY)
    
    @inbounds for ey in 1:NelY
        @inbounds for ex in 1:NelX
            eo = (ey-1)*NelX + ex

            # Properties of damage zone
            if (Domain_X*Domain-ThickX) >= ex*dxe >= ThickX && (dye <= ey*dye <= ThickY)   # ThickX: # ~distance from low boundary to fault zone low boundary
                mu .= rho2*vs2^2
                #  damage_idx[eo] = eo
            # Properties of host rock    
            else
                mu .= rho1*vs1^2
            end
            W[:,:,eo] = wgll2.*mu
        end
    end
    W
end

# Material properties for a trapezium damage zone
# Linear function: shape of trapezoid
# meaning of x,y: nodes of all elements
function line(x,y)
    P1 = [0 3e3]
    P2 = [-8e3 1.5e3]
    # f<0: host rock    f>0: damage rock  
    f = (y - P2[2]) - ((P1[2]-P2[2])/(P1[1]-P2[1]))*(x - P2[1])

    return f
end

# Set up trapezoidal rigidity
function rigid(x,y)
    # Rigidity: host rock and fault zone
    rho1::Float64 = 2670
    vs1::Float64 = 3464
    # inner fault zone
    rho2 = 0.6*rho1
    vs2 = 0.6*vs1
    # outer fault zone
    rho3 = 0.8*rho1
    vs3 = 0.8*vs1
    
    rhoglob::Array{Float64} = zeros(length(x))
    vsglob::Array{Float64} = zeros(length(x))

    for i in eachindex(x)
        if x[i] > -8e3
            if line(x[i],y[i]) < 0
                rhoglob[i] = rho3
                vsglob[i] = vs3
            else
                rhoglob[i] = rho1
                vsglob[i] = vs1
            end
        else
            rhoglob[i] = rho1
            vsglob[i] = vs1
        end

    end
    # properties of inner fault zone
    for i = 1:length(x)
        if y[i] < 0.25e3
            rhoglob[i] = rho2
            vsglob[i] = vs2
        end
    end

    # get the density and shear wave velocity of all elements
    return rhoglob, vsglob
end

function mat_trap(NelX, NelY, NGLL, iglob, M, dxe, dye, x,y, wgll2)
    # need to calculate the Mass matrix again
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    rho::Matrix{Float64} = zeros(NGLL, NGLL)
    W::Array{Float64,3} = zeros(NGLL, NGLL, NelX*NelY)
    rhoglob, vsglob = rigid(x,y)
    muglob = rhoglob.*(vsglob.^2)
    
    @inbounds for ey in 1:NelY
        @inbounds for ex in 1:NelX
            eo = (ey-1)*NelX + ex
            ig = iglob[:,:,eo]
            # heterogeneous mu and rho values
            mu[:,:] = muglob[ig]
            rho[:,:] = rhoglob[ig]
            
            W[:,:,eo] = wgll2.*mu
            M[ig] .+= wgll2.*rho*jac
        end
    end
    return M,W
end