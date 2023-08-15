include("$(@__DIR__)/compute_dot.jl")

x::Vector{Float64} = [1,2,3,4,5]
y::Vector{Float64} = [1,2,3,4,5]

compute_dot(x,y)