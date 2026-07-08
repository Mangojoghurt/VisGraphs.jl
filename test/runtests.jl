using VisGraphs
using Test

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

@testset "Basic functionality" begin
    x = generate_sine(10)

    hvg_edges = hvg(x)
    nvg_edges = nvg(x)

    @test length(hvg_edges) > 0
    @test length(nvg_edges) > 0
    @test hvg_edges != nvg_edges
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

@testset "Input validation" begin

    @test_throws ArgumentError hvg([1.0])
    @test_throws ArgumentError nvg([1.0])

    @test_throws ArgumentError hvg([1.0, Inf])
    @test_throws ArgumentError nvg([NaN, 1.0])

    @test_throws MethodError hvg([1.0, missing])
    @test_throws MethodError nvg([missing, 1.0])
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

@testset "Weighted graph input validation" begin

    @test_throws ArgumentError whvg([1.0])
    @test_throws ArgumentError wnvg([1.0])

    @test_throws ArgumentError whvg([1.0, Inf])
    @test_throws ArgumentError wnvg([NaN, 1.0])
end

@testset "Weighted plot functions" begin

    x = [1.0, 2.0, 1.5, 3.0]

    plt1 = plot_whvg(x)
    plt2 = plot_wnvg(x)

    @test plt1 !== nothing
    @test plt2 !== nothing
    @test typeof(plt1) == typeof(plt2)
end

@testset "adjacency_matrix — basic structure" begin

    x = [1.0, 2.0, 5.0]
    edges = nvg(x)
    A = adjacency_matrix(edges, length(x))

    @test size(A) == (3, 3)
    @test all(A[i, i] == 0 for i in 1:3)
    @test A == A'
    @test A[1, 2] == 1
    @test A[2, 3] == 1
    @test A[1, 3] == 1
    @test all(v in (0, 1) for v in A)
end

@testset "adjacency_matrix — HVG known example" begin

    x = [1.0, 2.0, 5.0]
    edges = hvg(x)
    A = adjacency_matrix(edges, length(x))

    @test A[1, 2] == 1
    @test A[2, 3] == 1
    @test A[1, 3] == 0
end

@testset "adjacency_matrix — works with weighted edges" begin

    x = [1.0, 2.0, 5.0]
    A = adjacency_matrix(wnvg(x), length(x))
    A_ref = adjacency_matrix(nvg(x), length(x))
    @test A == A_ref
end

@testset "adjacency_matrix — input validation" begin

    @test_throws ArgumentError adjacency_matrix([], 0)
    @test_throws ArgumentError adjacency_matrix([(1, 5)], 3)
end

@testset "degree_distribution — basic correctness" begin

    x = [1.0, 2.0, 5.0]

    edges = nvg(x)
    degrees, dist = degree_distribution(edges, length(x))

    @test length(degrees) == 3
    @test all(d == 2 for d in degrees)
    @test haskey(dist, 2)
    @test dist[2] ≈ 1.0

    edges_h = hvg(x)
    degrees_h, dist_h = degree_distribution(edges_h, length(x))

    @test degrees_h[1] == 1
    @test degrees_h[2] == 2
    @test degrees_h[3] == 1
    @test dist_h[1] ≈ 2/3
    @test dist_h[2] ≈ 1/3
end

@testset "degree_distribution — distribution sums to 1" begin

    x = generate_sine(30)

    for f in [hvg, nvg]
        edges = f(x)
        _, dist = degree_distribution(edges, length(x))
        @test sum(values(dist)) ≈ 1.0 atol=1e-12
    end
end

@testset "degree_distribution — degrees match adjacency matrix row sums" begin

    x = generate_noisy_sine(20)

    for f in [hvg, nvg]
        edges = f(x)
        n = length(x)
        degrees, _ = degree_distribution(edges, n)
        A = adjacency_matrix(edges, n)
        row_sums = [sum(A[i, :]) for i in 1:n]
        @test degrees == row_sums
    end
end

@testset "degree_distribution — input validation" begin

    @test_throws ArgumentError degree_distribution([], 0)
    @test_throws ArgumentError degree_distribution([(1, 5)], 3)
end

@testset "laplacian_matrix — basic properties" begin

    x = [1.0, 2.0, 5.0]
    edges = nvg(x)
    n = length(x)
    L = laplacian_matrix(edges, n)

    @test size(L) == (n, n)
    @test L == L'
    @test all(sum(L[i, :]) == 0 for i in 1:n)

    degrees, _ = degree_distribution(edges, n)
    @test all(L[i, i] == degrees[i] for i in 1:n)

    for i in 1:n, j in 1:n
        if i != j
            @test L[i, j] in (-1, 0)
        end
    end
end

@testset "laplacian_matrix — L = D - A" begin

    x = generate_sine(15)

    for f in [hvg, nvg]
        edges = f(x)
        n = length(x)
        A = adjacency_matrix(edges, n)
        D = zeros(Int, n, n)
        for i in 1:n
            D[i, i] = sum(A[i, :])
        end
        L = laplacian_matrix(edges, n)
        @test L == D - A
    end
end

@testset "laplacian_matrix — works with weighted edges" begin

    x = [1.0, 2.0, 5.0]
    L_weighted = laplacian_matrix(wnvg(x), length(x))
    L_unweighted = laplacian_matrix(nvg(x), length(x))
    @test L_weighted == L_unweighted
end

@testset "laplacian_matrix — input validation" begin

    @test_throws ArgumentError laplacian_matrix([], 0)
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

    using Random
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
