"""
Structure to handle the nonlinear complementarity problem in standard form:

   0 ≤ F(x) _|_ x ≥ 0
"""
struct NonLinearComplementarityProblem{T} <: AbstractComplementarityModel
  n # dimension of the problem
  F::Function
  Fx::AbstractVector{T}
  function NonLinearComplementarityProblem(n, F, Fx::AbstractVector{T}) where T
    return new{T}(n, F, Fx)
  end
end

Base.eltype(::NonLinearComplementarityProblem{T}) where {T} = T