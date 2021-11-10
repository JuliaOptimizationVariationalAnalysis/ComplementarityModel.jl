"""
Structure to handle the mixed complementarity problem in standard form:

   0 ≤ F(x) _|_ G(x) ≥ 0
"""
struct MixedComplementarityProblem{T} <: AbstractComplementarityModel
  n # dimension of the problem
  F::Function
  Fx::AbstractVector{T}
  G::Function
  Gx::AbstractVector{T}
  function MixedComplementarityProblem(
    n,
    F,
    Fx::AbstractVector{T},
    G,
    Gx::AbstractVector{T},
  ) where {T}
    return new{T}(n, F, Fx, G, Gx)
  end
end

Base.eltype(::MixedComplementarityProblem{T}) where {T} = T
