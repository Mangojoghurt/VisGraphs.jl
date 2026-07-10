using VisGraphs
using Test
using SparseArrays
using Random

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

@testset "adjacency_matrix — returns a sparse matrix" begin

    x = [1.0, 2.0, 5.0]
    A = adjacency_matrix(nvg(x), length(x))

    @test issparse(A)

    # sparsity change shouldn't affect any generic matrix behaviour
    @test size(A) == (3, 3)
    @test A == A'
    @test A[1, 2] == 1
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

@testset "degree_distribution — sorted keyword" begin

    x = [1.0, 2.0, 5.0]
    edges = hvg(x)  # degrees [1, 2, 1]

    degrees, dist = degree_distribution(edges, length(x); sorted=true)

    @test dist isa Vector
    @test issorted(first.(dist))
    @test dist[1] == (1 => 2/3)
    @test dist[2] == (2 => 1/3)
end

@testset "degree_distribution — default is still a Dict" begin

    x = [1.0, 2.0, 5.0]
    edges = hvg(x)

    degrees, dist = degree_distribution(edges, length(x))
    @test dist isa Dict{Int,Float64}
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

@testset "laplacian_matrix — returns a sparse matrix" begin

    x = [1.0, 2.0, 5.0]
    L = laplacian_matrix(nvg(x), length(x))

    @test issparse(L)
end

@testset "clustering_coefficient — fully connected triangle" begin

    # NVG on x=[1,2,5] connects all 3 nodes to each other
    x = [1.0, 2.0, 5.0]
    edges = nvg(x)

    local_c, global_c = clustering_coefficient(edges, length(x))

    @test all(c ≈ 1.0 for c in local_c)
    @test global_c ≈ 1.0
end

@testset "clustering_coefficient — path has zero clustering" begin

    # HVG on x=[1,2,5] is a path 1-2-3, no triangle
    x = [1.0, 2.0, 5.0]
    edges = hvg(x)

    local_c, global_c = clustering_coefficient(edges, length(x))

    @test all(c ≈ 0.0 for c in local_c)
    @test global_c ≈ 0.0
end

@testset "clustering_coefficient — input validation" begin

    @test_throws ArgumentError clustering_coefficient([], 0)
    @test_throws ArgumentError clustering_coefficient([(1, 5)], 3)
end

@testset "average_path_length — fully connected triangle" begin

    x = [1.0, 2.0, 5.0]
    edges = nvg(x)

    @test average_path_length(edges, length(x)) ≈ 1.0
end

@testset "average_path_length — path example" begin

    # HVG on x=[1,2,5] is a path 1-2-3: d(1,2)=1, d(2,3)=1, d(1,3)=2
    x = [1.0, 2.0, 5.0]
    edges = hvg(x)

    @test average_path_length(edges, length(x)) ≈ 4/3
end

@testset "average_path_length — input validation" begin

    @test_throws ArgumentError average_path_length([(1, 5)], 3)
    @test_throws ArgumentError average_path_length([], 1)
end

@testset "clustering_coefficient and average_path_length — random consistency" begin

    x = generate_noisy_sine(25)
    edges = nvg(x)
    n = length(x)

    local_c, global_c = clustering_coefficient(edges, n)

    @test length(local_c) == n
    @test all(0.0 ≤ c ≤ 1.0 for c in local_c)
    @test 0.0 ≤ global_c ≤ 1.0

    apl = average_path_length(edges, n)
    @test apl ≥ 1.0
    @test apl ≤ n - 1  # can't exceed the longest possible chain
end