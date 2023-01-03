using LogParadox
using Test
using Random

@testset "LogParadox.jl" begin
    # Write your tests here.
    @testset "gmam" begin
        Random.seed!(42)
        for i in 1:100
            xs = 1 .+ rand(1000)*i
            @test gm(xs) <= am(xs)
        end
    end

    @testset "mm" begin
        Random.seed!(42)
        for i in 1:100
            xs = 1 .+ rand(1000)*i
            ys = minmaxreplace(xs, i)
            @test gm(xs) <= gm(ys) <= am(ys) <= am(xs)
        end
    end
end
