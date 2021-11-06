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
  end
end