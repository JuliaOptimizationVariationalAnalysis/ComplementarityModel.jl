using ComplementarityModel
using Test

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
    @test (@allocated eval(fin)(Cx, Fx, Gx)) ≤ 32 # It is 32 for CFB
  end
end