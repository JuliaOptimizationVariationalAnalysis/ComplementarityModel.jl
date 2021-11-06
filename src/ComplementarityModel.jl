module ComplementarityModel

using LinearAlgebra, LinearOperators

include("Cfunctions.jl")

abstract type AbstractComplementarityModel end

include("lcp.jl")

export F, F!, G, G!, JF, JF!, JFv, JFv!, JG, JG!, JGv, JGv!
export LinearComplementarityProblem

end
