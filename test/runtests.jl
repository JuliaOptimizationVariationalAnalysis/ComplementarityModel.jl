using ComplementarityModel
using Random, Test
using LinearAlgebra, SparseArrays, LinearOperators

Random.seed!(1234)

@testset "ComplementarityModel.jl" begin
    # Write your tests here.
end

@testset "C-functions" begin
  n = 10
  for T in [Float16, Float64], field in [:Cmin, :CFB]
    Fx, Gx = rand(T, n), rand(T, n)
    @test eltype(eval(field)(Fx, Gx)) == T
    fin = Symbol(field, "!")
    Cx = similar(Fx)
    eval(fin)(Cx, Fx, Gx)
    @test (@allocated eval(fin)(Cx, Fx, Gx)) ≤ 32 # It is 32 for CFBf
    JFx, JGx = rand(T,n,n), rand(T,n,n)
    Cx = similar(JFx)
    Jfin = Symbol("J", field, "!")
    eval(Jfin)(Cx, Fx, Gx, JFx, JGx)
    @test (@allocated eval(Jfin)(Cx, Fx, Gx, JFx, JGx)) == 0
    Jfinv = Symbol("J", field, "v!")
    Jfv = Symbol("J", field, "v")
    w = similar(Fx)
    v = rand(T, n)
    eval(Jfv)(Fx, Gx, JFx, JGx, v)
    eval(Jfinv)(w, Fx, Gx, JFx, JGx, v)
    @test (@allocated eval(Jfinv)(w, Fx, Gx, JFx, JGx, v)) == 0
    @test w ≈ eval(Jfin)(Cx, Fx, Gx, JFx, JGx) * v atol = eps(T) * max(10, norm(v))
  end
end

@testset "LinearComplementarityProblem" begin
  T = Float64
  M = T[ 0 1 ; 1 0]
  Mlist = (
    M,
    sparse(M),
    LinearOperator(M),
  )
  q = rand(T, 2)
  for M in Mlist
    lcp = LinearComplementarityProblem(2, M, q)
    x = rand(T, 2)
    w = similar(x)
    @test F(lcp, x) == M * x + q
    F!(w, lcp, x)
    @test (@allocated F!(w, lcp, x)) ≤ 64 # not sure why though
    @test G(lcp, x) == x
    G!(w, lcp, x)
    @test (@allocated G!(w, lcp, x)) == 0
    @test JF(lcp, x) == M
    # @allocated JF(lcp, x)
    @test JG(lcp, x) == LinearAlgebra.I
    # @allocated JG(lcp, x)
    v = rand(T, 2)
    @test JFv(lcp, x, v) == M * v
    JFv!(w, lcp, x, v)
    @test (@allocated JFv!(w, lcp, x, v)) == 0
    @test JGv(lcp, x, v) == v
    JGv!(w, lcp, x, v)
    @test (@allocated JGv!(w, lcp, x, v)) == 0
  end
end
