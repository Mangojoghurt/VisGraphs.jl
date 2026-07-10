"""
    _degree_vector(edges, n::Integer)

Compute the degree of every node in `O(E)` time, where `E` is the number
of edges. Shared internal helper for `degree_distribution` and
`laplacian_matrix` so both only walk the edge list once.
"""
function _degree_vector(edges, n::Integer)

    n ≥ 1 || throw(ArgumentError("`n` must be at least 1."))

    degrees = zeros(Int, n)

    for edge in edges
        i, j = edge[1], edge[2]
        (1 ≤ i ≤ n && 1 ≤ j ≤ n) ||
            throw(ArgumentError("Edge ($i, $j) is out of range for n=$n."))
        degrees[i] += 1
        degrees[j] += 1
    end

    return degrees
end

"""
    _neighbor_lists(edges, n::Integer)

Build an adjacency list (one `Vector{Int}` of neighbors per node) from an
edge list, in `O(n + E)` time. Shared internal helper for
`clustering_coefficient` and `average_path_length`.
"""
function _neighbor_lists(edges, n::Integer)

    n ≥ 1 || throw(ArgumentError("`n` must be at least 1."))

    neighbors = [Int[] for _ in 1:n]

    for edge in edges
        i, j = edge[1], edge[2]
        (1 ≤ i ≤ n && 1 ≤ j ≤ n) ||
            throw(ArgumentError("Edge ($i, $j) is out of range for n=$n."))
        push!(neighbors[i], j)
        push!(neighbors[j], i)
    end

    return neighbors
end

"""
    adjacency_matrix(edges, n::Integer)

Construct an n×n **sparse** adjacency matrix from an edge list.

The graph is treated as undirected: each edge `(i, j)` contributes both
`A[i, j] = 1` and `A[j, i] = 1`. If edges are weighted tuples `(i, j, w)`,
the weight `w` is ignored.

Visibility graphs are typically sparse (each node connects to only a
small fraction of the others), so the matrix is built and returned as a
`SparseMatrixCSC`, using `O(E)` memory instead of the `O(n²)` a dense
matrix would need. All standard matrix operations (indexing, `size`,
`sum`, equality, etc.) work the same as on a dense matrix.

This function is intended for edge lists produced by `hvg`, `nvg`, `whvg`,
or `wnvg`.

# Example
```julia
x = [1.0, 2.0, 3.0, 1.5]
edges = nvg(x)
A = adjacency_matrix(edges, length(x))
```
"""
function adjacency_matrix(edges, n::Integer)

    n ≥ 1 || throw(ArgumentError("`n` must be at least 1."))

    is = Int[]
    js = Int[]

    for edge in edges
        i, j = edge[1], edge[2]
        (1 ≤ i ≤ n && 1 ≤ j ≤ n) ||
            throw(ArgumentError("Edge ($i, $j) is out of range for n=$n."))
        push!(is, i); push!(js, j)
        push!(is, j); push!(js, i)
    end

    vals = ones(Int, length(is))

    # combine keeps duplicate entries at 1 instead of summing them,
    # in case the edge list ever contains a repeated edge
    return sparse(is, js, vals, n, n, (a, b) -> one(a))
end

"""
    degree_distribution(edges, n::Integer; sorted::Bool=false)

Compute the degree sequence and degree distribution of a graph.

The degree of a node is defined as the number of incident edges.

# Keyword arguments
- `sorted::Bool=false`: if `true`, return the distribution as a
  `Vector{Pair{Int,Float64}}` sorted by ascending degree, instead of an
  unordered `Dict`. Useful when the result will be plotted or otherwise
  needs a deterministic order without an extra `sort(collect(dist))` step
  at the call site.

# Returns
- `degrees::Vector{Int}`: degree of each node.
- `distribution`: empirical degree distribution `P(k)`, as a
  `Dict{Int,Float64}` (default) or a sorted `Vector{Pair{Int,Float64}}`
  (when `sorted=true`).

# Example
```julia
x = generate_sine(50)
edges = nvg(x)

degrees, dist = degree_distribution(edges, length(x))
maximum(degrees)
sort(collect(dist))

# or, equivalently, ask for the sorted form directly:
degrees, dist_sorted = degree_distribution(edges, length(x); sorted=true)
```
"""
function degree_distribution(edges, n::Integer; sorted::Bool=false)

    degrees = _degree_vector(edges, n)

    # normalised frequency distribution P(k)
    distribution = Dict{Int,Float64}()
    for d in degrees
        distribution[d] = get(distribution, d, 0.0) + 1.0 / n
    end

    if sorted
        return degrees, Base.sort(collect(distribution); by=first)
    end

    return degrees, distribution
end

"""
    laplacian_matrix(edges, n::Integer)

Construct the combinatorial graph Laplacian `L = D - A`, where `D` is the
diagonal degree matrix and `A` is the (sparse) adjacency matrix.

The Laplacian is symmetric and positive semidefinite. The number of zero
eigenvalues corresponds to the number of connected components, and the
second-smallest eigenvalue (the Fiedler value) measures graph connectivity.

`D` is built as a `Diagonal` (storing only its `n` diagonal entries,
`O(n)`, rather than a dense `n×n` matrix), and combined with the sparse
`A`, so the result stays a sparse matrix — no `O(n²)` allocation is made
anywhere in this function.

# Example
```julia
x = generate_sine(20)
edges = nvg(x)

L = laplacian_matrix(edges, length(x))

using LinearAlgebra
eigvals(Symmetric(Matrix(L)))
```
"""
function laplacian_matrix(edges, n::Integer)

    A = adjacency_matrix(edges, n)
    degrees = _degree_vector(edges, n)
    D = Diagonal(degrees)

    return D - A
end

"""
    clustering_coefficient(edges, n::Integer)

Compute the local clustering coefficient of every node and the graph's
average (global) clustering coefficient.

For a node `v` with degree `k_v`, the local clustering coefficient is the
fraction of pairs of `v`'s neighbors that are themselves connected:

    C_v = (number of edges among neighbors of v) / (k_v choose 2)

Nodes with degree less than 2 have no well-defined ratio and are assigned
`C_v = 0.0`, following the convention of Watts & Strogatz (1998). The
global coefficient is the average of all `C_v`, taken over all `n` nodes.

This metric is discussed for visibility graphs in Lacasa et al. (2008).

# Returns
- `local_coeffs::Vector{Float64}`: local clustering coefficient of each node.
- `global_coeff::Float64`: average of `local_coeffs` over all `n` nodes.

# Example
```julia
x = generate_sine(30)
edges = nvg(x)
local_c, global_c = clustering_coefficient(edges, length(x))
```
"""
function clustering_coefficient(edges, n::Integer)

    neighbors = _neighbor_lists(edges, n)
    # sets give O(1) "is u a neighbor of v" checks in the triangle count below
    neighbor_sets = [Set(nb) for nb in neighbors]

    local_coeffs = zeros(Float64, n)

    for v in 1:n
        nb = neighbors[v]
        k = length(nb)
        k < 2 && continue

        links = 0
        for a in 1:k, b in a+1:k
            if nb[b] in neighbor_sets[nb[a]]
                links += 1
            end
        end

        local_coeffs[v] = 2 * links / (k * (k - 1))
    end

    global_coeff = sum(local_coeffs) / n

    return local_coeffs, global_coeff
end

"""
    average_path_length(edges, n::Integer)

Compute the average shortest-path length (in number of edges) between all
pairs of nodes, via breadth-first search from every node.

Visibility graphs built by `hvg`/`nvg` on a full time series are always
connected — every node is at least visible to its immediate neighbor —
so the average is well-defined over all `n*(n-1)/2` node pairs. If the
edge list passed in does describe a disconnected graph, an `ErrorException`
is thrown rather than silently ignoring unreachable pairs.

This metric is discussed for visibility graphs in Lacasa et al. (2008).

# Example
```julia
x = generate_sine(30)
edges = nvg(x)
L̄ = average_path_length(edges, length(x))
```
"""
function average_path_length(edges, n::Integer)

    n ≥ 2 || throw(ArgumentError("`n` must be at least 2 to define a path length."))

    neighbors = _neighbor_lists(edges, n)

    total = 0
    npairs = 0

    for s in 1:n
        dist = fill(-1, n)
        dist[s] = 0
        queue = [s]
        head = 1
        while head ≤ length(queue)
            u = queue[head]
            head += 1
            for v in neighbors[u]
                if dist[v] == -1
                    dist[v] = dist[u] + 1
                    push!(queue, v)
                end
            end
        end

        for t in s+1:n
            dist[t] == -1 &&
                error("Graph is disconnected; average path length is undefined between nodes $s and $t.")
            total += dist[t]
            npairs += 1
        end
    end

    return total / npairs
end
