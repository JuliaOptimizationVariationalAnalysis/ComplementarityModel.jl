module ComplementarityModel

using LinearAlgebra, LinearOperators, SparseArrays

include("Cfunctions.jl")

abstract type AbstractComplementarityModel end

include("lcp.jl")

"""
    LP_to_LCP(c, A, b)

Implements the self-embedding model for Linear Program

    min cᵀx s.t.  Ax=b, x ≥ 0   

from Roos, C., Terlaky, T., & Vial, J. P. (2005). Interior point methods for linear optimization.
"""
function LP_to_LCP(c::AbstractVector{T}, A, b::AbstractVector{T}) where {T}
  m, n = size(A)
  nn = n + m + 1
  M = spzeros(T, nn, nn)
  M[1:m, (m + 1):(m + n)] .= A
  M[1:m, (m + n + 1)] .= -b
  M[(m + 1):(m + n), 1:m] .= -A'
  M[(m + 1):(m + n), (m + n + 1)] .= c
  M[m + n + 1, 1:m] .= b
  M[m + n + 1, (m + 1):(m + n)] .= -c
  q = zeros(T, nn)
  return LinearComplementarityProblem(nn, M, q)
end

"""
    LP_to_feasible_LCP(c, A, b; μ = one(T))

Implements the self-embedding model for Linear Program

    min cᵀx s.t.  Ax=b, x ≥ 0   

from Roos, C., Terlaky, T., & Vial, J. P. (2005). Interior point methods for linear optimization.
The returned problem is such that ones(nn) * √μ is on the central path.
"""
function LP_to_feasible_LCP(c::AbstractVector{T}, A, b::AbstractVector{T}; μ = one(T)) where {T}
  n, m = size(A)
  nn = n + m + 2
  lcp = LP_to_LCP(c, A, b)
  Mbar = lcp.M
  r = ones(T, n + m + 1) - Mbar * ones(T, n + m + 1)
  M = [Mbar r; -r' 0]
  q = vcat(lcp.q, nn * √μ)
  return LinearComplementarityProblem(nn, M, q)
end

export LP_to_LCP, LP_to_feasible_LCP

include("ncp.jl")
include("mcp.jl")

export F, F!, G, G!, JF, JF!, JFv, JFv!, JG, JG!, JGv, JGv!
export LinearComplementarityProblem, NonLinearComplementarityProblem, MixedComplementarityProblem

end
