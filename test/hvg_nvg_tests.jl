using VisGraphs
using Test
using Random

@testset "Graph input validation" begin

    @test_throws ArgumentError hvg([1.0])
    @test_throws ArgumentError nvg([1.0])

    @test_throws ArgumentError hvg([1.0, Inf])
    @test_throws ArgumentError nvg([NaN, 1.0])

    @test_throws MethodError hvg([1.0, missing])
    @test_throws MethodError nvg([missing, 1.0])
end

@testset "Basic functionality" begin
    x = generate_sine(10)

    hvg_edges = hvg(x)
    nvg_edges = nvg(x)

    @test length(hvg_edges) > 0
    @test length(nvg_edges) > 0
    @test hvg_edges != nvg_edges
end

@testset "Large signal graph construction" begin
    x = generate_random(10_000)

    hvg_edges = hvg(x)
    nvg_edges = nvg(x)

    @test length(hvg_edges) > 0
    @test length(nvg_edges) > 0
end

@testset "HVG/NVG small examples" begin

    x = [1.0, 2.0, 5.0]

    # HVG
    hvg_edges = hvg(x)
    @test (1,2) in hvg_edges
    @test (2,3) in hvg_edges
    @test (1,3) ∉ hvg_edges

    # NVG
    nvg_edges = nvg(x)
    @test (1,2) in nvg_edges
    @test (2,3) in nvg_edges
    @test (1,3) in nvg_edges
end

@testset "Flat signals" begin

    x = fill(1.0, 5)

    vis_graphs = [hvg(x), nvg(x)]

    n = length(x)

    for edges in vis_graphs

        # flat signals should be chain-like
        @test length(edges) == n - 1

        for i in 1:n
            for j in i+1:n

                if j == i + 1
                    # adjacent edges must exist
                    @test (i, j) in edges
                else
                    # non-adjacent edges must NOT exist
                    @test (i, j) ∉ edges
                end
            end
        end
    end
end

@testset "Monotone signals" begin

    inc = collect(1:5)
    dec = collect(5:-1:1)

    expected_chain = [(i, i+1) for i in 1:4]

    for (name, x) in [("increasing", inc), ("decreasing", dec)]

        for (label, f) in [("HVG", hvg), ("NVG", nvg)]

            edges = f(x)

            # monotone signals should be chain-like
            @test length(edges) == length(x) - 1

            # adjacent edges must exist
            for e in expected_chain
                @test e in edges
            end
        end
    end
end

@testset "HVG vs NVG structural difference" begin

    x = [1.0, 3.0, 2.0, 4.0, 7.0]

    h = hvg(x)
    n = nvg(x)

    @test h != n
    @test length(h) < length(n)  # NVG usually denser
end

@testset "Plot functions" begin

    x = [1.0, 2.0, 1.5, 3.0]

    plt1 = plot_hvg(x)
    plt2 = plot_nvg(x)

    # ensure plot objects are returned
    @test plt1 !== nothing 
    @test plt2 !== nothing 

    # ensure correct type
    @test typeof(plt1) == typeof(plt2)
end

# ── Brute-force reference testing for the divide-and-conquer algorithms ────
#
# hvg()/nvg() use an O(n log n) divide-and-conquer implementation
# (src/core/internal.jl). These tests verify that implementation against
# a direct, unoptimized O(n^2)/O(n^3) reference that checks the visibility
# condition literally as stated in each function's docstring, across many
# random signals of different shapes.

"""
    _brute_hvg(x)

Reference HVG implementation: checks the horizontal visibility condition
`x[k] < min(x[i], x[j])` directly for every candidate pair `(i, j)`,
with no divide-and-conquer optimization.
"""
function _brute_hvg(x::AbstractVector{<:Real})
    n = length(x)
    edges = Tuple{Int,Int}[]
    for i in 1:n-1, j in i+1:n
        visible = all(x[k] < min(x[i], x[j]) for k in i+1:j-1)
        visible && push!(edges, (i, j))
    end
    return edges
end

"""
    _brute_nvg(x)

Reference NVG implementation: checks the natural visibility condition
`x[k] < x[i] + (x[j] - x[i]) * (k - i) / (j - i)` directly for every
candidate pair `(i, j)`, with no divide-and-conquer optimization.
"""
function _brute_nvg(x::AbstractVector{<:Real})
    n = length(x)
    edges = Tuple{Int,Int}[]
    for i in 1:n-1, j in i+1:n
        visible = all(x[k] < x[i] + (x[j] - x[i]) * (k - i) / (j - i) for k in i+1:j-1)
        visible && push!(edges, (i, j))
    end
    return edges
end

@testset "hvg/nvg — divide-and-conquer matches brute-force reference" begin

    rng = Random.MersenneTwister(42)

    signal_kinds = [
        n -> randn(rng, n),                                    # normal noise
        n -> rand(rng, n) .* 10 .- 5,                          # uniform noise
        n -> fill(3.14, n),                                    # constant
        n -> sort(rand(rng, n)),                               # monotone increasing
        n -> sort(rand(rng, n), rev=true),                     # monotone decreasing
        n -> (v = zeros(n); v[n ÷ 2 + 1] = 100.0; v),          # single spike
    ]

    for kind in signal_kinds, n in [2, 3, 5, 10, 20, 35]
        x = kind(n)
        @test hvg(x) == _brute_hvg(x)
        @test nvg(x) == _brute_nvg(x)
    end
end