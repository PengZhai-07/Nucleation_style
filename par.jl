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

function setParameters(FZdepth, halfwidth, res, T, alpha, multiple)

    LX::Int = 48e3  # depth dimension of rectangular domain
    LY::Int = 32e3 # off fault dimenstion of rectangular domain

    NelX::Int = 30*res # no. of elements in x
    NelY::Int = 20*res # no. of elements in y

    dxe::Float64 = LX/NelX   #	Size of one element along X
    print(dxe)
    dye::Float64 = LY/NelY   #	Size of one element along Y
    print(dye)
    Nel::Int = NelX*NelY     # Total no. of elements

    P::Int = 4		#	Lagrange polynomial degree
    NGLL::Int = P + 1 #	No. of Gauss-Legendre-Lobatto nodes
    FltNglob::Int = NelX*(NGLL - 1) + 1  # No. of GLL nodes in X direction

    # Jacobian for global -> local coordinate conversion
    # For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
    # is linear, its jacobian is constant
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    coefint1::Float64 = jac/dx_dxi^2         # dye/dxe
    coefint2::Float64 = jac/dy_deta^2         # dxe/dye
    print(coefint1)
    print(coefint2)

    #..................
    # TIME PARAMETERS
    #..................

    yr2sec::Int = 365*24*60*60

    Total_time::Int = T*yr2sec     # Set the total time for simulation here

    CFL::Float64 = 0.6	#	Courant stability number     c*(dt/dx) <= 1

    IDstate::Int = 2    #   State variable equation type: aging law

    # Some other variables used to define the time steps: see koneko's function
    dtincf::Float64 = 1.2      # useless
    gamma_::Float64 = pi/4
    dtmax::Int = 400 * 24 * 60*60		# Is this useful to determine the time step

    #...................
    # MEDIUM PROPERTIES
    #...................

    # default: host rock!!
    rho1::Float64 = 2670
    vs1::Float64 = 3464
    
    # # The entire medium has low rigidity
    # rho1::Float64 = 2670
    # vs1::Float64 = 0.6*3464

    # the initial property of fualt damage zone: fault zone evolution!!!
    rho2::Float64 = 2670
    vs2::Float64 = sqrt(alpha)*vs1

    mu = rho1*vs1^2
    print("The shear modulus of shostrock is",mu)     # the default value is about 32GPa

    # without viscosity damping
    ETA = 0.

    # Low velocity layer dimensions
    ThickX::Float64 = LX - ceil(FZdepth/dxe)*dxe   # ~FZdepth m deep
    ThickY::Float64 = ceil(halfwidth/dye)*dye     # ~ default halfwidth value: 0.25 km wide
    # when the resolution is low, the halfwidth of fault damage zone can not be two small 

    #.......................
    # EARTHQUAKE PARAMETERS
    #.......................

    Vpl::Float64 = 35e-3/yr2sec	#	Plate loading rate   unit: m/seconds
    # frictional parameters along the fault line (X direction)
    fo::Vector{Float64} = repeat([0.6], FltNglob)       #	Reference friction coefficient
    Vo::Vector{Float64} = repeat([1e-6], FltNglob)		#	Reference velocity 'Vo'  unit: m/s
    xLf::Vector{Float64} = repeat([0.008], FltNglob)    #	Dc (Lc) = 8 mm

    Vthres::Float64 = 0.001     # unit: m/s  if slip rate is higher than this value, earthquake happens
    Vevne::Float64 = Vthres     # redefine the velocity threshold!!

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
    # gg = iglob[:,:,1];    # 2D matrix
    # println(gg)
    x = x .- LX     # +x direction is upward (-48km ~ 0km)
    # println("size of x=", size(x))
    #println("size of iglob=", size(iglob))
    nglob::Int = length(x)    # no. of all GLL nodes  

    # The derivatives of the Lagrange Polynomials were pre-tabulated
    # xgll = location of the GLL nodes inside the reference segment [-1,1] : same to XGLL in MeshBox.jl
    # Using the GetGLL function
    xgll::Vector{Float64}, wgll::Vector{Float64}, H::Matrix{Float64} = GetGLL(NGLL)
    wgll2::SMatrix{NGLL,NGLL,Float64} = wgll*wgll'       # define a sparse Matrix from FEMSparse
    # println(wgll2)


    #.............................
    #   OUTPUT RECEIVER LOCATIONS
    #.............................
    # For now, it saves slip, sliprate, and stress at the nearest node specified.
    # My coordinates are weird, might change them later.
    # x coordinate = along dip fault length (always -ve(value) below the free surface)   upward is positive
    # y coordinate = off-fault distance (+ve)     right is positive              

    ## we could output the theoretical seismogram on this points: only SH waves
    ## and SH waves will only affect the shear stress on the fault line 
    x_out = [6.0, 6.0, 6.0, 6.0, 6.0, 6.0].*(-1e3)  # x coordinate of receiver
    y_out = [66.0, 130.0, 198.0, 250.0, 330.0, 396.0]     # y coordinate of receiver   # 
    #  n_receiver = length(x_receiver) # number of receivers
    x_out, y_out, out_seis, dist = FindNearestNode(x_out, y_out, x, y)

    #.................
    # Initialization
    #.................

    # For internal forces: Local contributions to the stiffness matrix K
    W::Array{Float64,3} = zeros(NGLL, NGLL, Nel)
    
    # Global Mass Matrix
    M::Vector{Float64} = zeros(nglob)
    # println("Dimension of M=", size(M))

    # Mass + Damping matrix
    #  MC::Vector{Float64} = zeros(nglob)

    # Assemble mass and stiffness matrix
    # Mass Assembly: M = wgll2.*rho*jac

    a , M, dt::Float64, muMax = Massemble!(NGLL, NelX, NelY, dxe, dye,
        ThickX, ThickY, rho1, vs1, rho2, vs2, iglob, M, x, y, jac)
    
    # println(M[1:100])
    # println(a)          # total number of elements
    # Material properties for a narrow rectangular damaged zone of half-thickness ThickY and depth ThickX: shear modulus
    W = MaterialProperties(NelX, NelY, NGLL, dxe, dye, ThickX, ThickY, wgll2, rho1, vs1, rho2, vs2)
    # W[:,:,eo] = wgll2.*mu

    # Material properties for trapezoid damaged zone: need to calculate mass matrix again!! may not calculate the dt again,
    # because the velocity of host rock(vs1: the bigger value) decides the timestep!!!(without change of mesh!!!) 
    # M, W =  mat_trap(NelX, NelY, NGLL, iglob, M, dxe, dye, x, y, wgll2)

    # Stiffness Assembly: compute and output the whole large stiffness matrix: about (154401*154401)
    Ksparse::SparseMatrixCSC{Float64} = Kassemble(NGLL, NelX, NelY, dxe,dye, nglob, iglob, W)
    #print(size(Ksparse))

    # Damage Indexed Kdam
    # fault damage zone evolution
    did = damage_indx!(ThickX, ThickY, dxe, dye, NGLL, NelX, NelY, iglob)
    #println("index of GLL nodes in fault damage zone: ", did) 

    #  return Ksparse, Kdam, iglob
    #  Kdam[Kdam .> 1.0] .= 1.0

    # Time solver variables
    # current dafault dt meets c*(dt/dx) <= 1 by function 'Massemble'
    dt = CFL*dt
    dtmin = dt
    half_dt = 0.5*dtmin
    # half_dt_sq = 0.5*dtmin^2

    #......................
    # Boundary conditions :  L->B  R->T  T->R  B->L
    #......................

    # Bottom boundary: absorbing
    BcBC::Vector{Float64}, iBcB::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, 
                        rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'B') 
    
    # Top Boundary = free surface: nothing to do
    BcTC::Vector{Float64}, iBcT::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY,
                        rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'T')

    # Right Boundary: absorbing
    BcRC::Vector{Float64}, iBcR::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, 
                        rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'R')

    # Mass matrix at boundaries: absorbing
    #  Mq = M[:]
    M[iBcB] .= M[iBcB] .+ half_dt*BcBC
    M[iBcR] .= M[iBcR] .+ half_dt*BcRC
    #  M[iBcT] .= M[iBcT] .+ half_dt*BcRT

    # for boundary conditions on fault line   : 'L' 
    # iFlt: index of GLL nodes on the fault!!
    FltL::Vector{Float64}, iFlt::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, 
                       rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'L')

    # what is FltZ?  FltZ = jac1D(dy_deta)*wgll*rho1/dt_min  used in NRsearch
    # M = wgll2.*rho1 (damage zone).*jac (dx_dxi*dy_deta)
    # FltL = dx_dxi.*wgll.*1 (impedance)
    # half_dt = 0.5*dtmin
    FltZ::Vector{Float64} = M[iFlt]./FltL/half_dt * 0.5   #  times 0.5 due to the symmetry 
    print(FltZ)
    # X (vertical) of all GLL nodes at fault surface
    FltX::Vector{Float64} = x[iFlt]   
    #println("# X(vertical) of all GLL nodes at dynamic fault surface: ", FltX)    
    
    #......................
    # Initial Conditions
    #......................
    cca::Vector{Float64}, ccb::Vector{Float64}, a_b = fricDepth(FltX)   # rate-state friction parameters
    
    Seff::Vector{Float64} = SeffDepth(FltX, multiple)       # default effective normal stress: 10MPa

    Snormal::Vector{Float64} = SnormalDepth(FltX)       # effective normal stress
    # print(Snormal)

    SSpp::Vector{Float64} = SSppDepth(FltX)       # effective normal stress
    # print(SSpp)

    tauo::Vector{Float64} = tauDepth(FltX, multiple)        # initial shear stress

    # Kelvin-Voigt Viscosity : one technical method to increase the convergence rate
    Nel_ETA::Int = 0   # not used! 
    if ETA != 0
        Nel_ETA = NelX
        x1 = 0.5*(1 .+ xgll')
        eta_taper = exp.(-pi*x1.^2)
        eta = ETA*dt*repeat([eta_taper], NGLL)
    else
        Nel_ETA = 0
    end

    # Compute XiLF(largest slip in one timestep!!) used in timestep calculation: constrained by friction law!
    # quasi-static scheme
    XiLf::Vector{Float64} = XiLfFunc!(LX, FltNglob, gamma_, xLf, muMax, cca, ccb, Seff)

    # Find nodes that do not belong to the fault (off-fault GLL nodes )
    FltNI::Vector{Int} = deleteat!(collect(1:nglob), iFlt)

    # Compute diagonal of K
    #  diagKnew::Vector{Float64} = KdiagFunc!(FltNglob, NelY, NGLL, Nel, coefint1, coefint2, iglob, W, H, Ht, FltNI)

    # Fault boundary: global indices where fault within 24 km: boundary between dynamic fault and creeping fault  
    fbc = reshape(iglob[:,1,:], length(iglob[:,1,:]))   #convert the index of all left(fault) boundary GLL nodes in all elements into 1-D vector
    # println("fbc=", fbc[1:10])
    # println(findall(x .== -24e3)[1])    # the point on the fault at the depth of 24km
    idx = findall(fbc .== findall(x .== -24e3)[1] - 1)[1]
    #println("idx=", idx)
    FltIglobBC::Vector{Int} = fbc[1:idx]     # GLL nodes within creeping fault (>24 km)  with repeated nodes
    
    # Display important parameters
    # println("Total number of GLL nodes on fault: ", FltNglob)
    println("Total number of GLL nodes on fault: ", length(iFlt))
    println("Average node spacing on fault: ", LX/(FltNglob-1), " m")
    println("Effective normal stress distribution on fault: \n ", Seff, " MPa")
    println("friction parameter a-b on fault:\n", a_b)
    println("ThickX: ", ThickX, " m")
    println("ThickY: ", ThickY, " m")
    @printf("dt: %1.09f s\n", dt)   # minimal timestep during coseismic stage

    return params_int(Nel, FltNglob, yr2sec, Total_time, IDstate, nglob),
            params_float(ETA, Vpl, Vthres, Vevne, dt, mu, ThickY),
            # arrary = vector
            params_farray(fo, Vo, xLf, M, BcBC, BcRC, FltL, FltZ, FltX, cca, ccb, Seff, Snormal, SSpp, 
            tauo, XiLf, x_out, y_out),
            params_iarray(iFlt, iBcB, iBcR, FltIglobBC, FltNI, out_seis), 
            Ksparse, iglob, NGLL, wgll2, nglob, did

end



# the sequence of above parameters should be the same with the following
struct params_int{T<:Int}
    # Domain size
    Nel::T          # total number of elements in the  2D model
    FltNglob::T     # total number of GLL nodes on the fault line  

    # Time parameters
    yr2sec::T       # how many seconds in a year
    Total_time::T   # Total simulation time
    IDstate::T      # the type of the friction law

    # Fault setup parameters
    nglob::T        # total number of GLL nodes in the 2D model

end

struct params_float{T<:AbstractFloat}
    # Jacobian for global -> local coordinate conversion
    #  jac::T
    #  coefint1::T
    #  coefint2::T
    # shear modulus
     
    ETA::T    # do not use the 

    # Earthquake parameters
    Vpl::T     # plate motion rate used in this model
    Vthres::T   # velocity threshold for coseismic phase
    Vevne::T    # same meaning with Vthres

    # Setup parameters
    dt::T        # timestep based on CFL creterion   
    mu::T        # shear modulus of host rock
    ThickY::T     # real halfwidth of damage zone
end

struct params_farray{T<:Vector{Float64}}
    fo::T         # reference friction coefficient on the whole fault 
    Vo::T         # reference velocity on the whole fault
    xLf::T        # characteristic slip distance Dc

    M::T      # mass vector of the whole model

    BcBC::T    # boundary condition on the Bottom(mass)
    BcRC::T    # boundary condition on the right(mass)

    FltL::T    # boundary condition on the fault(mass)
    FltZ::T     # 
    FltX::T   # real depth of all nodes on the fault 

    cca::T    # a of RSF 
    ccb::T    # b of RSF 
    Seff::T   # effective normal stress
    Snormal::T   # initial normal stress
    SSpp::T   # initial steady state pore pressure
    tauo::T    # intial shear stress
 
    XiLf::T    # maximum slip in a timestep to constrain the length of timestep (based on friction law)
    #  diagKnew::T

    xout::T   # output the seismogram
    yout::T   # output the seismogram
end

struct params_iarray{T<:Vector{Int}}
    iFlt::T        # index of GLL nodes on the fault
    iBcB::T        # index of GLL nodes on the bottom boundary
    iBcR::T        # index of GLL nodes on the right boundary
    FltIglobBC::T  # index of GLL nodes within the creep fault (twice of the shared nodes)
    FltNI::T       # index of off-fault GLL nodes
    out_seis::T    # index of off-fault GLL nodes which are nearest to the predefined output locations
end

# required timestepping for low slip velocities(quasi-static)
# Calculate XiLf used in computing the final timestep
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

        if Xithf*Xith[j] > Ximax
            XiLf[j] = Ximax*xLf[j]          # 0.5*DC
        else
            XiLf[j] = Xithf*Xith[j]*xLf[j]
        end
    end
    
    return XiLf
end
