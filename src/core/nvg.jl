"""
    nvg(x)

Construct a Natural Visibility Graph (NVG) from a time series `x`.

Two nodes `i < j` are connected if every intermediate point lies below
the straight line segment joining them, i.e.:

    x[k] < x[i] + (x[j] - x[i]) * (k - i) / (j - i)
    for all i < k < j

The resulting graph encodes convex visibility relations in the time series.

Assumptions:
- samples are indexed on an evenly spaced time grid (t = 1, 2, ..., n)
- the index position is used as the time coordinate

For irregularly sampled data, users must resample or provide explicit
timestamps (not currently supported).

The implementation uses a divide-and-conquer algorithm, based on the approach used in the Python package `ts2vg`.

Returns a sorted edge list `(i, j)` with `i < j`.
"""
function nvg(x::AbstractVector{<:Real})
    _validate(x)
    edges = Tuple{Int,Int}[]
    sizehint!(edges, 6 * length(x))
    _nvg_core!(edges, x, 1, length(x))
    return sort!(edges)
end

"""
    wnvg(x)

Construct a Weighted Natural Visibility Graph (WNVG) from a time series `x`.

Extends the NVG by assigning each edge a geometric weight equal to the
angle of visibility between two points:

    w(i, j) = atan(x[j] - x[i], j - i)

The weight encodes both vertical and horizontal separation between nodes.

Assumptions:
- samples are indexed on an evenly spaced time grid (t = 1, 2, ..., n)
- the index position is used as the time coordinate

For irregularly sampled data, users must resample or provide explicit
timestamps (not currently supported).

Returns a vector of weighted edges `(i, j, w)`.
"""
function wnvg(x::AbstractVector{<:Real})
    base  = nvg(x) 
    W     = float(eltype(x))
    edges = Tuple{Int,Int,W}[]
    sizehint!(edges, length(base))
    
    @inbounds for (i, j) in base
        weight = atan(W(x[j] - x[i]), W(j - i))
        push!(edges, (i, j, weight))
    end
    return edges
end

"""
    plot_nvg(x)

Plot a time series together with its Natural Visibility Graph (NVG).

The time series is shown as a line plot, and NVG edges are drawn as
straight-line connections between mutually visible points.

Assumptions:
- samples are indexed on an evenly spaced time grid (t = 1, 2, ..., n)
- the index position is used as the time coordinate

For irregularly sampled data, users must resample or provide explicit
timestamps (not currently supported).

Returns a `Plots.Plot` object.
"""
function plot_nvg(x::AbstractVector{<:Real})
    edges = nvg(x)

    T = float(eltype(x))
    xs = T[]
    ys = T[]

    for (i, j) in edges
        push!(xs, i, j, NaN)
        push!(ys, x[i], x[j], NaN)
    end

    plt = plot(xs, ys;
        color=:gray,
        alpha=0.5,
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
    plot_wnvg(x)

Plot a time series together with its Weighted Natural Visibility Graph (WNVG).

Edges are colored according to their normalized angular weight across the graph.

Assumptions:
- samples are indexed on an evenly spaced time grid (t = 1, 2, ..., n)
- the index position is used as the time coordinate

For irregularly sampled data, users must resample or provide explicit
timestamps (not currently supported).

Returns a `Plots.Plot` object.
"""
function plot_wnvg(x::AbstractVector{<:Real})
    edges = wnvg(x)

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