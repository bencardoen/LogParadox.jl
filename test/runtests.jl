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

    @testset "rep" begin
        Random.seed!(42)
        replacers = [rep_min!, rep_max!, rep_minmax!, rep_rand!]
        for i in 1:100
            xs = 1 .+ rand(1000)*i
            for replacer in replacers
                ys = transform_steps_replace(xs, 100, replacer)
                @test am(ys) != am(xs)
                @test gm(ys) != gm(xs)
            end
        end
    end


    @testset "sm" begin
        Random.seed!(42)

        xs = 1 .+ rand(1000)*10

        s = smooth(xs, 5)
        t = smooth(xs, 10)
        @test ! all(s .== xs[1:length(s)])
        @test ! all(t .== s[1:length(t)])
    end

    @testset "repadd" begin
        Random.seed!(42)
        for i in 1:100
            xs = 1 .+ rand(1000)*i
            ys = transform_steps_replace(xs, 100)
            @test am(ys) < am(xs)
            @test gm(ys) > gm(xs)
        end
    end
end
