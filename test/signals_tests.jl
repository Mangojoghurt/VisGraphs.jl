using VisGraphs
using Test
using Random

@testset "Signal generators" begin

    s = generate_sine(20)
    r = generate_random(20)
    ns = generate_noisy_sine(20)

    @test length(s) == 20
    @test length(r) == 20
    @test length(ns) == 20

    @test all(isfinite, s)
    @test all(isfinite, r)
    @test all(isfinite, ns)
end

@testset "Signal generator length validation" begin

    for n in (-10, 0, 1)
        @test_throws ArgumentError generate_sine(n)
        @test_throws ArgumentError generate_random(n)
        @test_throws ArgumentError generate_noisy_sine(n)
    end
end

@testset "Signal generator RNG reproducibility" begin

    rng1 = MersenneTwister(42)
    rng2 = MersenneTwister(42)

    r1 = generate_random(20; rng=rng1)
    r2 = generate_random(20; rng=rng2)

    @test r1 == r2

    rng1 = MersenneTwister(42)
    rng2 = MersenneTwister(42)

    ns1 = generate_noisy_sine(20; rng=rng1)
    ns2 = generate_noisy_sine(20; rng=rng2)

    @test ns1 == ns2
end

@testset "Signal generator RNG independence" begin

    r1 = generate_random(20; rng=MersenneTwister(42))
    r2 = generate_random(20; rng=MersenneTwister(43))

    @test r1 != r2

    ns1 = generate_noisy_sine(20; rng=MersenneTwister(42))
    ns2 = generate_noisy_sine(20; rng=MersenneTwister(43))

    @test ns1 != ns2
end