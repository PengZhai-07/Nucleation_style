using Base
import Random: seed!
using SparseArrays
using LinearAlgebra
using AlgebraicMultigrid
import AlgebraicMultigrid as AMG
using BenchmarkTools

# Poisson matrix also makes the point, but the second-level A with Ruge-Stuben coarsening is a bit too dense,
# so the gap between sparse and dense solves is not significant enough.
# Here use a very sparse random matrix that lead to a sparse-enough second level.
# A = poisson((64, 64))
seed!(0)
A = sprand(8000, 8000, 1e-4) * 0.1 + I  #  add to diagonal to avoid singular matrix
b = A * ones(size(A)[1])

ruge_stuben(A, max_levels=10)  # fast if coarsest A is very small
@time ml = ruge_stuben(A, max_levels=2)  # takes 5 seconds, mostly spends on `inv(Matrix(ml.final_A))`
@btime AMG._solve(ml, b, maxiter=10);  # takes 16 ms


import SuiteSparse.UMFPACK: UmfpackLU
import AlgebraicMultigrid: CoarseSolver

# note: UMFPACK only supports Float64 or ComplexF64
# not sure how to best restrict the types
struct Splu <: CoarseSolver     # an abstract supertype
    LU::UmfpackLU
    Splu(A) = new(lu(A))  # here A is only a new parameter
end
Base.show(io::IO, p::Splu) = print(io, "  Splu")    # translation p into a display function

function (p::Splu)(x,b)  
    x[:] = p.LU \ b
end



# println(dump(Splu))
ml_sp = ruge_stuben(A, max_levels=2, coarse_solver=Splu)  
@btime ruge_stuben(A, max_levels=2, coarse_solver=Splu)  # takes 6 ms, ~1000x faster than before!
@btime AMG._solve(ml_sp, b, maxiter=10);  # takes 9 ms, 2x faster than before