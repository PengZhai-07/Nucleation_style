function compute_dot(DX::Vector{Float64}, DY::Vector{Float64})
    @assert length(DX) == length(DY)     # prints a diagnostic message to stderr and aborts the program if expression is false
    n = length(DX)
    incx = incy = 1
    product = ccall((:ddot_, "libLAPACK"),
            Float64,
            (Ref{Int32}, Ptr{Float64},Ref{Int32}, Ptr{Float64}, Ref{Int32}),
            n, DX, incx, DY, incy)
    return product
end
