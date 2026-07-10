using VisGraphs
using Test

@testset "Weighted graph input validation" begin

    @test_throws ArgumentError whvg([1.0])
    @test_throws ArgumentError wnvg([1.0])

    @test_throws ArgumentError whvg([1.0, Inf])
    @test_throws ArgumentError wnvg([NaN, 1.0])

    @test_throws MethodError whvg([1.0, missing])
    @test_throws MethodError wnvg([missing, 1.0])
end

@testset "Weighted graphs return correct type" begin

    x = [1.0, 2.0, 5.0]

    wh = whvg(x)
    wn = wnvg(x)

    @test eltype(wh) == Tuple{Int,Int,Float64}
    @test eltype(wn) == Tuple{Int,Int,Float64}
end

@testset "Weighted graphs have same edges as unweighted" begin

    x = generate_sine(20)

    # strip weights and compare edge sets to unweighted versions
    @test [(i,j) for (i,j,_) in whvg(x)] == hvg(x)
    @test [(i,j) for (i,j,_) in wnvg(x)] == nvg(x)
end

@testset "Weighted graph weights are finite" begin

    x = generate_noisy_sine(30)

    for (_,_,w) in whvg(x)
        @test isfinite(w)
    end

    for (_,_,w) in wnvg(x)
        @test isfinite(w)
    end
end

@testset "Weight formula is correct" begin

    # for two points (1, x1) and (2, x2) the angle is atan(x2-x1, 1)
    x = [0.0, 1.0]

    wh = whvg(x)
    wn = wnvg(x)

    @test length(wh) == 1
    @test length(wn) == 1

    expected_weight = atan(1.0, 1.0)  # π/4

    @test wh[1][3] ≈ expected_weight
    @test wn[1][3] ≈ expected_weight
end

@testset "Weighted plot functions" begin

    x = [1.0, 2.0, 1.5, 3.0]

    plt1 = plot_whvg(x)
    plt2 = plot_wnvg(x)

    @test plt1 !== nothing
    @test plt2 !== nothing
    @test typeof(plt1) == typeof(plt2)
end