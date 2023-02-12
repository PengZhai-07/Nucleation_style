# Linear interpolation function
function Int1D(P1, P2, val)	
	Line = P1[1] .+ ( (P2[1] - P1[1])/((P2[2] - P1[2])).*(val .- P1[2]) )	
	return Line
end

# define the rate and state friction parameter in all 48km long fault 
# Compute rate-state friciton with depth
function fricDepth(FltX, asp_a, asp_b, matrix_a, Domain, multiple_matrix, multiple_asp)
    
    FltNglob = length(FltX)    # number of GLL nodes on fault
    
    asp_ab = asp_a-asp_b      # a-b value in asperity
    matrix_ab = matrix_a - asp_b      # a-b value in background matrix

    ccb::Array{Float64} = repeat(asp_b, FltNglob)    # b is always a constant
    a_b::Array{Float64} = repeat(matrix_ab, FltNglob)     # tempory a_b equals matrix_ab

    # setup the transiton for kinematic fault and RSF fault
    # [a-b, depth]   key points of friction coefficient change
    fP1 = [0.024, 0e3]  
    fP2 = [0.024, -40e3*Domain/12]
    fP3 = [matrix_ab, -40e3*Domain/6]
    fP4 = [matrix_ab, -40e3*Domain*5/6]
    fP5 = [0.024, -40e3*Domain*11/12]
    fP6 = [0.024, -40e3*Domain]

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

    # setup the distribution of asperities and background matrix 
    L_fault = fP3[2]-fP4[2]      # length of fault(asperity + background matrix)
    N::Int = 2e6
    n::Int = 4     # 1 asperity every n matrix
    cell_size = L_fault/N     # cell size is about 62.5m
    N_group::Int = floor(N/(n+1))
    N_remain = N - N_group*(n+1)    # put it at the beginning of the group

    cca .= a_b .+ asp_b      # so a is variable and b is a constant in all depth

    NS = multiple_matrix*10e6    #  tempory Seff equals matrix
    Seff::Array{Float64} = repeat([NS], FltNglob)

    for i = 1:N_group
        index_depth = findall(abs(fP3[2])+ N_remain*cell_size+(i-1)*(n+1)*cell_size .<= abs.(FltX) .<= abs(fP3[2])+ N_remain*cell_size+(i-1)*(n+1)*cell_size+ cell_size)
        a_b[index_depth] .= asp_ab
        Seff[index_depth] .= multiple_asp*10e6
    end
    tauo::Array{Float64} = Seff.*0.6

    return cca, ccb, a_b, Seff, tauo
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
