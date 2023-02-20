# Linear interpolation function
function Int1D(P1, P2, val)	
	Line = P1[1] .+ ( (P2[1] - P1[1])/((P2[2] - P1[2])).*(val .- P1[2]) )	
	return Line
end

# define the rate and state friction parameter in all 48km long fault 
# Compute rate-state friciton with depth
function fricDepth(FltX, asp_a, asp_b, matrix_a, Domain, multiple_matrix, multiple_asp, matrix_asp_ratio,N,asperity_number)
    
    FltNglob = length(FltX)    # number of GLL nodes on fault
    
    asp_ab = asp_a-asp_b      # a-b value in asperity
    matrix_ab = matrix_a - asp_b      # a-b value in background matrix

    cca::Array{Float64} = repeat([asp_a], FltNglob)    # tempory a equals matrix_ab
    ccb::Array{Float64} = repeat([asp_b], FltNglob)    # b is always a constant
    a_b::Array{Float64} = repeat([matrix_ab], FltNglob)     # tempory a_b equals matrix_ab

    # setup the transiton for kinematic fault and RSF fault
    # [a-b, depth]   key points of friction coefficient change
    fP1 = [0.047, 0e3]  
    fP2 = [0.047, -Domain_X*Domain/8]
    fP3 = [matrix_ab, -Domain_X*Domain/4]
    fP4 = [matrix_ab, -Domain_X*Domain*3/4]
    fP5 = [0.047, -Domain_X*Domain*7/8]
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

    # setup the distribution of asperities and background matrix 
    L_fault = fP3[2]-fP4[2]      # length of fault(asperity + background matrix)
    cell_size = L_fault/N     # cell size should be about 100 m
    N_remain::Int = N - asperity_number*(matrix_asp_ratio+1)   # put it at the beginning of the group
    
    println("Cell size: ", cell_size)
    println("Total number of groups: ", N_group)
    println("Remain cells: ", N_remain)

    NS = multiple_matrix*10e6    #  tempory Seff equals matrix
    Seff::Array{Float64} = repeat([NS], FltNglob)

    for i = 1:N_group

        index_depth = findall(abs(fP3[2])+ N_remain*cell_size+(i-1)*(matrix_asp_ratio+1)*cell_size .<= abs.(FltX) .<= abs(fP3[2])+ N_remain*cell_size+(i-1)*(matrix_asp_ratio+1)*cell_size+cell_size)
        println("asperity GLL node index:", index_depth)
        # set the transiton
        a_b[index_depth] .= asp_ab

        index_transition_left = collect((index_depth[1] - floor(Int, length(index_depth)/2)):index_depth[1])
        println("The transition points: ",index_transition_left)
        fP_a = [matrix_ab, FltX[index_transition_left[1]]]
        fP_b = [asp_ab, FltX[index_transition_left[end]]]
        a_b[index_transition_left] .= Int1D(fP_a, fP_b, FltX[index_transition_left])

        index_transition_right = collect(index_depth[end]:(index_depth[end] + floor(Int, length(index_depth)/2)))
        println("The transition points: ",index_transition_right)
        fP_a = [asp_ab, FltX[index_transition_right[1]]]
        fP_b = [matrix_ab, FltX[index_transition_right[end]]]
        a_b[index_transition_right] .= Int1D(fP_a, fP_b, FltX[index_transition_right])
       
        Seff[index_depth] .= multiple_asp*10e6
        
    end
    tauo::Array{Float64} = Seff.*0.6

    cca .= a_b .+ asp_b      # so a is variable and b is a constant in all depth

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
