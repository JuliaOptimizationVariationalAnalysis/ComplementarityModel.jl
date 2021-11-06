module ComplementarityModel

using LinearAlgebra, LinearOperators

include("Cfunctions.jl")

abstract type AbstractComplementarityModel end

include("lcp.jl")
include("ncp.jl")
include("mcp.jl")

export F, F!, G, G!, JF, JF!, JFv, JFv!, JG, JG!, JGv, JGv!
export LinearComplementarityProblem, NonLinearComplementarityProblem, MixedComplementarityProblem

end
