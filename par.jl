#######################################################################
#	PARAMETER FILE: SET THE PHYSICAL PARAMETERS FOR THE SIMULATION
#######################################################################
include("$(@__DIR__)/src/GetGLL.jl")		 #	Polynomial interpolation
include("$(@__DIR__)/src/MeshBox.jl")		 # 	Build 2D mesh
include("$(@__DIR__)/src/FindNearestNode.jl")   #	Nearest node for output
include("$(@__DIR__)/src/Massemble.jl")       #   Assemble mass matrix
include("$(@__DIR__)/src/Kassemble.jl")      #   Assemble stiffness matrix
include("$(@__DIR__)/src/MaterialProperties.jl")		 # 	Build 2D mesh
#  include("$(@__DIR__)/trapezoidFZ/Assemble.jl") #   Gaussian fault zone assemble
include("$(@__DIR__)/src/damageEvol.jl")   #    Stiffness index of damaged medium
include("$(@__DIR__)/src/BoundaryMatrix.jl")    #	Boundary matrices
include("$(@__DIR__)/src/initialConditions/defaultInitialConditions.jl")


function setParameters(FZdepth, halfwidth, res, T)

    LX::Int = 48e3  # depth dimension of rectangular domain
    LY::Int = 30e3 # off fault dimenstion of rectangular domain

    NelX::Int = 30*res # no. of elements in x
    NelY::Int = 20*res # no. of elements in y

    dxe::Float64 = LX/NelX   #	Size of one element along X
    dye::Float64 = LY/NelY   #	Size of one element along Y
    Nel::Int = NelX*NelY     # Total no. of elements

    P::Int = 4		#	Lagrange polynomial degree
    NGLL::Int = P + 1 #	No. of Gauss-Legendre-Lobatto nodes
    FltNglob::Int = NelX*(NGLL - 1) + 1  # No. of GLL nodes in X direction

    # Jacobian for global -> local coordinate conversion
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    coefint1::Float64 = jac/dx_dxi^2          # dye/dxe
    coefint2::Float64 = jac/dy_deta^2         # dxe/dye

    #..................
    # TIME PARAMETERS
    #..................

    yr2sec::Int = 365*24*60*60

    Total_time::Int = T*yr2sec     # Set the total time for simulation here

    CFL::Float64 = 0.6	#	Courant stability number     c*(dt/dx) <= 1

    IDstate::Int = 2    #   State variable equation type: aging law

    # Some other time variables used in the loop
    dtincf::Float64 = 1.2
    gamma_::Float64 = pi/4
    dtmax::Int = 400 * 24 * 60*60		# 400 days   


    #...................
    # MEDIUM PROPERTIES
    #...................

    # default
    rho1::Float64 = 2670
    vs1::Float64 = 3464

    # The entire medium has low rigidity
    #  rho1::Float64 = 2500
    #  vs1::Float64 = 0.6*3464

    # elastic parameters of damage zone?
    rho2::Float64 = 2670
    vs2::Float64 = 1.00*vs1

    # without viscosity
    ETA = 0.

    # Low velocity layer dimensions
    ThickX::Float64 = LX - ceil(FZdepth/dxe)*dxe   # ~FZdepth m deep
    ThickY::Float64 = ceil(halfwidth/dye)*dye     # ~ default halfwidth value: 0.25 km wide
    # when the resolution is low, the halfwidth of fault damage zone can not be two small 

    #.......................
    # EARTHQUAKE PARAMETERS
    #.......................

    Vpl::Float64 = 35e-3/yr2sec	#	Plate loading rate   unit: m/seconds
    # frictional parameters along the fault (X direction)
    fo::Vector{Float64} = repeat([0.6], FltNglob)       #	Reference friction coefficient
    Vo::Vector{Float64} = repeat([1e-6], FltNglob)		#	Reference velocity 'Vo'  unit: m/s
    xLf::Vector{Float64} = repeat([0.008], FltNglob)    #	Dc (Lc) = 8 mm

    Vthres::Float64 = 0.005     # unit: m/s  if slip rate is higher than this value, earthquake happens
    Vevne::Float64 = Vthres     #??

    #-----------#
    #-----------#
    # SETUP
    #-----------#
    #-----------#

    #....................
    # 2D Mesh generation
    #....................
    # global node index of the (i,j) th GLL node internal to the e-th element.
    iglob::Array{Int,3}, x::Vector{Float64}, y::Vector{Float64} =
                        MeshBox!(NGLL, Nel, NelX, NelY, FltNglob, dxe, dye)
    x = x .- LX     # +x direction is upward (-48km ~ 0km)
    nglob::Int = length(x)    # no. of all GLL nodes  

    # The derivatives of the Lagrange Polynomials were pre-tabulated
    # xgll = location of the GLL nodes inside the reference segment [-1,1]   same to XGLL in MeshBox.jl
    # Using the GetGLL function
    xgll::Vector{Float64}, wgll::Vector{Float64}, H::Matrix{Float64} = GetGLL(NGLL)
    wgll2::SMatrix{NGLL,NGLL,Float64} = wgll*wgll'       # define a sparse Matrix from FEMSparse
    
    #.............................
    #   OUTPUT RECEIVER LOCATIONS
    #.............................
    # For now, it saves slip, sliprate, and stress at the nearest node specified.
    # My coordinates are weird, might change them later.
    # x coordinate = along dip fault length (always -ve(value) below the free surface)   upward is positive
    # y coordinate = off-fault distance (+ve)     right is positive              

    x_out = [6.0, 6.0, 6.0, 6.0, 6.0, 6.0].*(-1e3)  # x coordinate of receiver
    y_out = [66.0, 130.0, 198.0, 250.0, 330.0, 396.0]     # y coordinate of receiver   # 
    #  n_receiver = length(x_receiver) # number of receivers

    x_out, y_out, out_seis, dist = FindNearestNode(x_out, y_out, x, y)

    #.................
    # Initialization
    #.................

    # For internal forces
    # W::Array{Float64,3} = zeros(NGLL, NGLL, Nel)

    # Global Mass Matrix
    M::Vector{Float64} = zeros(nglob)

    # Mass + Damping matrix
    #  MC::Vector{Float64} = zeros(nglob)

    # Assemble mass and stiffness matrix
    # Mass Assembly
    M, dt::Float64, muMax = Massemble!(NGLL, NelX, NelY, dxe, dye,
        ThickX, ThickY, rho1, vs1, rho2, vs2, iglob, M, x, y, jac)

    # Material properties for a narrow rectangular damaged zone of
    # half-thickness ThickY and depth ThickX
    W = MaterialProperties(NelX, NelY, NGLL, dxe, dye, 
        ThickX, ThickY, wgll2, rho1, vs1, rho2, vs2)

    # Material properties for trapezoid damaged zone
    #  M, W =  mat_trap(NelX, NelY, NGLL, iglob, M, dxe, dye, x, y, wgll2)

    # Stiffness Assembly
    Ksparse::SparseMatrixCSC{Float64} = Kassemble(NGLL, NelX, NelY, dxe,dye, nglob, iglob, W)

    # Damage Indexed Kdam
    # fault damage zone evolution
    did = damage_indx!(ThickX, ThickY, dxe, dye, NGLL, NelX, NelY, iglob)
    #println("index of GLL nodes in fault damage zone: ", did) 

    #  return Ksparse, Kdam, iglob
    #  Kdam[Kdam .> 1.0] .= 1.0

    # Time solver variables
    # current dafault dt meets c*(dt/dx) <= 1  from function 'Massemble'
    dt = CFL*dt
    dtmin = dt
    half_dt = 0.5*dtmin
    half_dt_sq = 0.5*dtmin^2

    #......................
    # Boundary conditions :  L->B  R->T  T->R  B->L
    #......................

    # Bottom boundary
    BcBC::Vector{Float64}, iBcB::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, 
                        rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'B') 
    
    # Top Boundary = free surface: nothing to do
    #  BcTC, iBcT = BoundaryMatrix(P, wgll, iglob, 'T')

    # Right Boundary
    BcRC::Vector{Float64}, iBcR::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, 
                        rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'R')

    # Mass matrix at boundaries
    #  Mq = M[:]
    M[iBcB] .= M[iBcB] .+ half_dt*BcBC
    M[iBcR] .= M[iBcR] .+ half_dt*BcRC
    #  M[iBcT] .= M[iBcT] .+ half_dt*BcRT

    # for initial conditions on fault line   : 'L' 
    FltL::Vector{Float64}, iFlt::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, 
                       rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'L')
    # iFlt: index of GLL nodes on the fault!!
    FltZ::Vector{Float64} = M[iFlt]./FltL /half_dt * 0.5   # specific meanings?
    # X(vertical) of all GLL nodes at fault surface
    FltX::Vector{Float64} = x[iFlt]   
    #println("# X(vertical) of all GLL nodes at dynamic fault surface: ", FltX)    
    
    #......................
    # Initial Conditions
    #......................
    cca::Vector{Float64}, ccb::Vector{Float64} = fricDepth(FltX)   # rate-state friction parameters
    Seff::Vector{Float64} = SeffDepth(FltX)       # effective normal stress
    tauo::Vector{Float64} = tauDepth(FltX)        # initial shear stress

    # Kelvin-Voigt Viscosity : one kind of initial condition?
    Nel_ETA::Int = 0
    if ETA !=0
        Nel_ETA = NelX
        x1 = 0.5*(1 .+ xgll')
        eta_taper = exp.(-pi*x1.^2)
        eta = ETA*dt*repeat([eta_taper], NGLL)

    else
        Nel_ETA = 0
    end

    # Compute XiLF(largest slip in one timestep!!) used in timestep calculation
    XiLf::Vector{Float64} = XiLfFunc!(LX, FltNglob, gamma_, xLf, muMax, cca, ccb, Seff)

    # Find nodes that do not belong to the fault (off-fault GLL nodes )
    FltNI::Vector{Int} = deleteat!(collect(1:nglob), iFlt)

    # Compute diagonal of K
    #  diagKnew::Vector{Float64} = KdiagFunc!(FltNglob, NelY, NGLL, Nel, coefint1, coefint2, iglob, W, H, Ht, FltNI)

    # Fault boundary: indices where fault within 24 km 
    fbc = reshape(iglob[:,1,:], length(iglob[:,1,:]))    # convert the index of all left boundasy GLL nodes in all elements  into 1-D vector
    idx = findall(fbc .== findall(x .== -24e3)[1] - 1)[1]
    FltIglobBC::Vector{Int} = fbc[1:idx]
    println("Total number of GLL nodes on fault boundary: ", length(FltIglobBC)) # ?
    
    # Display important parameters
    println("Total number of GLL nodes on fault: ", FltNglob)
    println("Total number of GLL nodes on fault: ", length(iFlt))
    println("Average node spacing: ", LX/(FltNglob-1), " m")
    println("ThickX: ", ThickX, " m")
    println("ThickY: ", ThickY, " m")
    @printf("dt: %1.09f s\n", dt)   # minimal timestep during coseismic stage

    return params_int(Nel, FltNglob, yr2sec, Total_time, IDstate, nglob),
            params_float(ETA, Vpl, Vthres, Vevne, dt),
            params_farray(fo, Vo, xLf, M, BcBC, BcRC, FltL, FltZ, FltX, cca, ccb, Seff, tauo, XiLf, x_out, y_out),
            params_iarray(iFlt, iBcB, iBcR, FltIglobBC, FltNI, out_seis), 
            Ksparse, iglob, NGLL, wgll2, nglob, did

end

struct params_int{T<:Int}
    # Domain size
    Nel::T
    FltNglob::T

    # Time parameters
    yr2sec::T
    Total_time::T
    IDstate::T

    # Fault setup parameters
    nglob::T

end

struct params_float{T<:AbstractFloat}
    # Jacobian for global -> local coordinate conversion
    #  jac::T
    #  coefint1::T
    #  coefint2::T

    ETA::T

    # Earthquake parameters
    Vpl::T
    Vthres::T
    Vevne::T

    # Setup parameters
    dt::T
end

struct params_farray{T<:Vector{Float64}}
    fo::T
    Vo::T
    xLf::T

    M::T

    BcBC::T
    BcRC::T

    FltL::T
    FltZ::T
    FltX::T

    cca::T
    ccb::T
    Seff::T
    tauo::T

    XiLf::T
    #  diagKnew::T

    xout::T
    yout::T
end
struct params_iarray{T<:Vector{Int}}
    iFlt::T
    iBcB::T
    iBcR::T
    FltIglobBC::T
    FltNI::T
    out_seis::T
end

# Calculate XiLf used in computing the timestep
function XiLfFunc!(LX, FltNglob, gamma_, xLf, muMax, cca, ccb, Seff)

    hcell = LX/(FltNglob-1)   # average interval of GLL nodes on fault
    Ximax = 0.5
    Xithf = 1

    Xith:: Vector{Float64} = zeros(FltNglob)
    XiLf::Vector{Float64} = zeros(FltNglob)

    #  @inbounds for j = 1:FltNglob
    @inbounds for j = 1:FltNglob

        # Compute time restricting parameters
        expr1 = -(cca[j] - ccb[j])/cca[j]
        expr2 = gamma_*muMax/hcell*xLf[j]/(cca[j]*Seff[j])  

        ro = expr2 - expr1

        if (0.25*ro*ro - expr2) >= 0
            Xith[j] = 1/ro
        else
            Xith[j] = 1 - expr1/expr2
        end

        # For each node, compute slip that node cannot exceed in one timestep
        # slip during one time step could not exceed 0.5*DC, so that we can restrict the length of the timestep!!

        if Xithf*Xith[j] > Ximax
            XiLf[j] = Ximax*xLf[j]
        else
            XiLf[j] = Xithf*Xith[j]*xLf[j]
        end
    end

    return XiLf
end


