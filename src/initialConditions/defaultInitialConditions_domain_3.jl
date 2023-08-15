# Linear interpolation function
function Int1D(P1, P2, val)	
	Line = P1[1] .+ ( (P2[1] - P1[1])/((P2[2] - P1[2])).*(val .- P1[2]) )	
	return Line
end

# define the rate and state friction parameter in all 48km long fault 
# Compute rate-state friciton with depth
function fricDepth(FltX, Domain)
    
    FltNglob = length(FltX)    # number of GLL nodes on fault
    
    # Friction with depth
    cca::Array{Float64} = repeat([0.015], FltNglob)
    ccb::Array{Float64} = repeat([0.019], FltNglob)    # b is always a constant

    a_b = cca - ccb     # -0.004 is the initial value of a-B
    # [a-b, depth]   key points of friction coefficient change
    fP1 = [0.047, 0e3]   # fP1 = [-0.003, 0e3]
    fP2 = [0.024, -17.5e3+0.25e3]
    fP3 = [-0.004, -17.5e3]
    fP4 = [-0.004, -22.5e3]
    fP5 = [0.024, -22.5e3-0.25e3]
    fP6 = [0.047, -Domain_X*Domain]

    # Return a vector I of the indices or keys of A
    fric_depth1 = findall(abs.(FltX) .<= abs(fP2[2]))
    fric_depth2 = findall(abs(fP2[2]) .< abs.(FltX) .<= abs(fP3[2]))
    fric_depth3 = findall(abs(fP3[2]) .< abs.(FltX) .<= abs(fP4[2]))
    fric_depth4 = findall(abs(fP4[2]) .< abs.(FltX) .<= abs(fP5[2]))
    fric_depth5 = findall(abs(fP5[2]) .< abs.(FltX) .<= abs(fP6[2]))

    a_b[fric_depth1] .= Int1D(fP1, fP2, FltX[fric_depth1])
    a_b[fric_depth2] .= Int1D(fP2, fP3, FltX[fric_depth2])
    a_b[fric_depth3] .= Int1D(fP3, fP4, FltX[fric_depth3])
    a_b[fric_depth4] .= Int1D(fP4, fP5, FltX[fric_depth4])
    a_b[fric_depth5] .= Int1D(fP5, fP6, FltX[fric_depth5])

    #  cca[fric_depth4] .= Int1D(fP4, fP5, FltX[fric_depth4]) .+ 0.0001
    cca .= ccb .+ a_b      # so a is variable and b is a constant in all depth
    #  ccb .= cca .- a_b

    return cca, ccb, a_b

end

# Effective normal stress
function SeffDepth(FltX, multiple)

    FltNglob = length(FltX)
    NS = multiple*10e6
    Seff::Array{Float64} = repeat([NS], FltNglob)
    # sP1 = [0.2*NS 0]
    # sP2 = [NS -2e3]         # constant normal stress below 2 km
    # Seff_depth = findall(abs.(FltX) .<= abs(sP2[2]))
    # Seff[Seff_depth] = Int1D(sP1, sP2, FltX[Seff_depth])

    return Seff

end

# Shear stress
function tauDepth(FltX, multiple)

    FltNglob = length(FltX)
    NS = multiple*10e6
    tauo::Array{Float64} = repeat([0.6*NS], FltNglob)
    tP1 = [0.45*NS 0]      
    #tP1 = [0.01e6 0]    
    tP2 = [0.45*NS -17.5e3 + 0.25e3]
    #  tP2 = [30e6 -0.5e3]
    tP3 = [0.6*NS -17.5e3]
    tP4 = [0.6*NS -22.5e3]
    tP5 = [0.45*NS -22.5e3 - 0.25e3]
    tP6 = [0.45*NS -Domain_X*Domain]

    tau_depth1 = findall(abs.(FltX).<=  abs(tP2[2]))
    tau_depth2 = findall(abs(tP2[2]) .< abs.(FltX) .<= abs(tP3[2]))
    tau_depth3 = findall(abs(tP3[2]) .< abs.(FltX) .<= abs(tP4[2]))
    tau_depth4 = findall(abs(tP4[2]) .< abs.(FltX) .<= abs(tP5[2]))
    tau_depth5 = findall(abs(tP5[2]) .< abs.(FltX) .<= abs(tP6[2]))

    tauo[tau_depth1] .= 0.45*NS
    tauo[tau_depth2] = Int1D(tP2, tP3, FltX[tau_depth2])
    tauo[tau_depth4] = Int1D(tP4, tP5, FltX[tau_depth4])
    tauo[tau_depth5] .= 0.45*NS

    return tauo
end














# test functions
# Initial normal stress: linear dependent
function SnormalDepth(FltX)

    FltNglob = length(FltX)
    Snormal::Array{Float64} = repeat([0], FltNglob)
    sP1 = [10e6 0]
    sP2 = [970e6 -48e3]  
    Snormal_depth = findall(abs.(FltX) .<= abs(sP2[2]))
    Snormal[Snormal_depth] = Int1D(sP1, sP2, FltX[Snormal_depth])

    return Snormal

end

# Initial pore pressure: linear dependent
function SSppDepth(FltX)

    FltNglob = length(FltX)
    SSpp::Array{Float64} = repeat([0], FltNglob)
    sP1 = [10e6 0]
    sP2 = [50e6 -5e3]  
    sP3 = [920e6 -48e3]
    SSpp_depth1 = findall(abs.(FltX) .<= abs(sP3[2]))
    SSpp_depth2 = findall(abs.(FltX) .<= abs(sP2[2]))
    SSpp[SSpp_depth1] = Int1D(sP2, sP3, FltX[SSpp_depth1])
    SSpp[SSpp_depth2] = Int1D(sP1, sP2, FltX[SSpp_depth2])

    return SSpp

end
