"""
Structure to handle the linear complementarity problem in standard form:

   0 ≤ Mx + q _|_ x ≥ 0
"""
struct LinearComplementarityProblem{T, S} <: AbstractComplementarityModel
  n # dimension of the problem
  M::S # matrix, sparsematrix or operator of size (n, n)
  q::AbstractVector{T} # vector of size n
  function LinearComplementarityProblem(n, M, q::AbstractVector{T}) where T
    return new{T, typeof(M)}(n, M, q)
  end
end

function F(lcp::LinearComplementarityProblem, x)
  Fx = similar(lcp.q)
  return F!(Fx, lcp, x)
end

function F!(Fx, lcp::LinearComplementarityProblem, x)
  mul!(Fx, lcp.M, x)
  Fx .+= lcp.q
  return Fx
end

function JF(lcp::LinearComplementarityProblem, x)
  JFx = similar(lcp.M)
  return JF!(JFx, lcp, x)
end

function JF(lcp::LinearComplementarityProblem{T, S}, x) where {T, S <: LinearOperator}
  JFx = lcp.M
  return JF!(JFx, lcp, x)
end

function JF!(JFx, lcp::LinearComplementarityProblem, x)
  JFx .= lcp.M
  return JFx
end

function JF!(JFx, lcp::LinearComplementarityProblem{T, S}, x) where {T, S <: LinearOperator}
  return lcp.M # broadcast is not implemented for LinearOperators
end

function JFv(lcp::LinearComplementarityProblem, x, v)
  w = similar(x)
  return JFv!(w, lcp, x, v)
end

function JFv!(w, lcp::LinearComplementarityProblem, x, v)
  mul!(w, lcp.M, v)
  return w
end

function G(lcp::LinearComplementarityProblem, x)
  Gx = similar(lcp.q)
  return G!(Gx, lcp, x)
end

function G!(Gx, lcp::LinearComplementarityProblem, x)
  Gx .= x
  return Gx
end

function JG(lcp::LinearComplementarityProblem, x)
  JGx = similar(lcp.M)
  return JG!(JGx, lcp, x)
end

function JG(lcp::LinearComplementarityProblem{T, S}, x) where {T, S <: LinearOperator}
  JGx = I # similar is not implemented for LinearOperators
  return JG!(JGx, lcp, x)
end

function JG!(JGx, lcp::LinearComplementarityProblem, x)
  JGx = I
  return JGx
end

function JGv(lcp::LinearComplementarityProblem, x, v)
  w = similar(x)
  return JGv!(w, lcp, x, v)
end

function JGv!(w, lcp::LinearComplementarityProblem, x, v)
  w .= v
  return w
end
