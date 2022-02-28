# Function to find the mesh node that is closest to the requested 'nseis' locations
function FindNearestNode(xin, yin, X, Y)
    # xin:  x coordinate of receiver
    # yin:  y coordinate of receiver
    nseis = length(xin)
    dist = zeros(nseis)
    iglob::Vector{Int} = zeros(nseis)

    for k = 1:nseis
        iglob[k] = argmin( (X .- xin[k]).^2 + (Y .- yin[k]).^2 )    # find the index of k nearest points
        dist[k] = minimum( (X .- xin[k]).^2 + (Y .- yin[k]).^2 )    # find the distance of k nearest points
    end

    dist = sqrt.(dist)
    xout = X[iglob]
    yout = Y[iglob]

    return xout, yout, iglob, dist

end
