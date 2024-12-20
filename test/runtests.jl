using LogParadox
using Test
using Random
using StatsBase
using DataFrames

@testset "LogParadox.jl" begin
    # Write your tests here.
    @testset "gmam" begin
        Random.seed!(42)
        for i in 1:100
            xs = 1 .+ rand(1000)*i
            @test gm(xs) <= am(xs)
        end
    end

    @testset "gemi" begin
        Random.seed!(42)
        X, Y = 512, 512
        ## Sizes of objects
        sizes = [1, 3, 9, 27]
        ## Frequencies per celltype
        freq_a = [300, 100, 30, 7]
        freq_b = [240, 147, 30, 4]
        Na = sum(freq_a)
        Nb = sum(freq_b)
        pa = freq_a ./ Na
        pb = freq_b ./ Nb
        S = 100
        ac, bc, ima, imb, ga, gb, GSa, GSb = generate_images_from_markov_chains(pa, pb, sizes, sizes; X=X, Y=Y, S=S, matchstate=4)
        @test length(ac) == length(bc)
        @test length(bc) == 4
        @test ima != imb
        freq_a = [300, 100, 30, 7]
        freq_b = [240, 147, 30, 0]
        Na = sum(freq_a)
        Nb = sum(freq_b)
        pa = freq_a ./ Na
        pb = freq_b ./ Nb
        S = 100
        ac, bc, ima, imb, ga, gb, GSa, GSb = generate_images_from_markov_chains(pa, pb, sizes, sizes; X=X, Y=Y, S=S, matchstate=4)
        @test length(ac) == length(bc)
        @test length(bc) == 4
        @test length(GSa) == length(GSb)
        @test ima != imb
    end

    @testset "mm" begin
        Random.seed!(42)
        for i in 1:100
            xs = 1 .+ rand(1000)*i
            ys = minmaxreplace(xs, i)
            @test gm(xs) <= gm(ys) <= am(ys) <= am(xs)
        end
    end

    @testset "checkdf" begin 
        Random.seed!(2)
        df = DataFrame(rand(50, 10), :auto)
        df.label = df.x9 .> .05
        rdf = check_dataframe(df)
        @test size(rdf) == (10, 8)
        @test any(rdf[:, "paradox?"])
    end

    @testset "testlp" begin
        X = [2, 4, 6, 13.0]
        Y = [5.5]
        ys = vcat(X, Y)
        # am(X)
        # gm(X)
        # am(ys)
        # gm(ys)
        gx, gy, ax, ay, p = check_paradox(X, ys)
        @test p
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

    @testset "repadd" begin
        Random.seed!(42)
        for i in 1:100
            xs = 1 .+ rand(1000)*i
            ys = transform_steps(xs, 100)
            @test am(ys) != am(xs)
            @test gm(ys) != gm(xs)
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


    @testset "scenario" begin

        X, Y = 512, 512

        ## Sizes of objects
        sizes = [1, 3, 9, 27]
        ## Frequencies per celltype
        freq_a = [300, 100, 30, 7]
        freq_b = [240, 147, 30, 4]

        Na = sum(freq_a)
        Nb = sum(freq_b)

        pa = freq_a ./ Na
        pb = freq_b ./ Nb

        S = 100
        Random.seed!(42)
        a_counts = nothing
        b_counts = nothing
        while true
            rs = rand(S)
            entries_a = to_entries(rs, pa)
            entries_b = to_entries(rs, pb)
            a_counts = countmap(entries_a)
            b_counts = countmap(entries_b)
            if length(a_counts)==4 && length(b_counts)==4
                if a_counts[4] == b_counts[4]
                    break
                end
            end
        end
        imga = generate_image(a_counts, sizes, X, Y)
        imgb = generate_image(b_counts, sizes, X, Y)
        @test imga != imgb

    end
end
