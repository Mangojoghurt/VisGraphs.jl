"""
    hvg(x)

Construct a Horizontal Visibility Graph (HVG) from a time series `x`.

Two nodes `i < j` are connected if all intermediate values satisfy
the horizontal visibility condition:

    x[k] < min(x[i], x[j])  for all i < k < j

The resulting graph encodes the structure of visibility relationships
in the time series.

Assumptions:
- samples are indexed on an evenly spaced time grid (t = 1, 2, ..., n)
- the index position is used as the time coordinate

For irregularly sampled data, users must resample or provide explicit
timestamps (not currently supported).

The implementation uses a divide-and-conquer algorithm with worst-case O(N^2)
complexity, based on the approach used in the Python package `ts2vg`.

Returns a sorted edge list `(i, j)` with `i < j`.
"""
function hvg(x::AbstractVector{<:Real})
    _validate(x)
    edges = Tuple{Int,Int}[]
    # Pre-allocate to minimize dynamic memory resizing (JUML H2 optimization)
    sizehint!(edges, 4 * length(x))
    _hvg_core!(edges, x, 1, length(x))
    return sort!(edges)
end

"""
    whvg(x)

Construct a Weighted Horizontal Visibility Graph (WHVG) from a time series `x`.

Extends the HVG by assigning each edge a geometric weight equal to the
angle of visibility between two connected points:

    w(i, j) = atan(x[j] - x[i], j - i)

The weight encodes both amplitude difference and temporal separation.

Assumptions:
- samples are indexed on an evenly spaced time grid (t = 1, 2, ..., n)
- the index position is used as the time coordinate

For irregularly sampled data, users must resample or provide explicit
timestamps (not currently supported).

Returns a vector of weighted edges `(i, j, w)`.
"""
function whvg(x::AbstractVector{<:Real})
    base  = hvg(x)
    W     = float(eltype(x)) # Respects Float32/Float64 inputs natively
    edges = Tuple{Int,Int,W}[]
    sizehint!(edges, length(base))
    
    @inbounds for (i, j) in base
        weight = atan(W(x[j] - x[i]), W(j - i))
        push!(edges, (i, j, weight))
    end
    return edges
end

"""
    plot_hvg(x)

Plot a time series together with its Horizontal Visibility Graph (HVG).

The time series is shown as a line plot, while edges are drawn as
semi-transparent connections between visible nodes.

Assumptions:
- samples are indexed on an evenly spaced time grid (t = 1, 2, ..., n)
- the index position is used as the time coordinate

For irregularly sampled data, users must resample or provide explicit
timestamps (not currently supported).

Returns a `Plots.Plot` object.
"""
function plot_hvg(x::AbstractVector{<:Real})
    edges = hvg(x)

    T = float(eltype(x))
    xs = T[]
    ys = T[]

    for (i, j) in edges
        push!(xs, i, j, NaN)
        push!(ys, x[i], x[j], NaN)
    end

    plt = plot(xs, ys;
        color=:gray,
        alpha=0.6,
        label=false,
        xlabel="t",
        ylabel="x(t)",
        title="Horizontal Visibility Graph"
    )

    plot!(plt, 1:length(x), x;
        lw=2,
        label="time series",
        color=:blue
    )

    return plt
end

"""
    plot_whvg(x)

Plot a time series together with its Weighted Horizontal Visibility Graph (WHVG).

Edges are colored according to their angular weight, normalized across
all edges.

Assumptions:
- samples are indexed on an evenly spaced time grid (t = 1, 2, ..., n)
- the index position is used as the time coordinate

For irregularly sampled data, users must resample or provide explicit
timestamps (not currently supported).

Returns a `Plots.Plot` object.
"""
function plot_whvg(x::AbstractVector{<:Real})
    edges = whvg(x)

    T = float(eltype(x))

    weights = T[w for (_, _, w) in edges]
    wmin, wmax = extrema(weights)
    wrange = wmax - wmin

    xs = T[]
    ys = T[]
    colors = RGB{T}[]

    for (i, j, w) in edges
        c = wrange ≈ zero(T) ? T(0.5) : (w - wmin) / wrange

        push!(xs, i, j, NaN)
        push!(ys, x[i], x[j], NaN)

        color = RGB{T}(c, zero(T), one(T) - c)
        push!(colors, color, color, color)
    end

    plt = plot(xs, ys;
        color=colors,
        alpha=0.6,
        label=false,
        xlabel="t",
        ylabel="x(t)",
        title="Weighted Horizontal Visibility Graph"
    )

    plot!(plt, 1:length(x), x;
        lw=2,
        label="time series",
        color=:blue
    )

    return plt
end