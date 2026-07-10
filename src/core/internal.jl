"""
    _argmax_range(x, left, right)

Return the index of the maximum value in `x[left:right]`.

If multiple entries share the maximum value, return the leftmost index.

This is an internal helper used in divide-and-conquer visibility graph
construction.
"""
@inline function _argmax_range(x::AbstractVector{<:Real}, left::Int, right::Int)
    idx  = left
    best = @inbounds x[left]
    @inbounds for k in left+1:right
        xk = x[k]
        if xk > best
            best = xk
            idx  = k
        end
    end
    return idx
end

"""
    _hvg_core!(edges, x, left, right)

Construct Horizontal Visibility Graph (HVG) edges for the subrange
`x[left:right]` and append them to `edges`.

The algorithm uses a divide-and-conquer strategy with an explicit stack
to avoid recursion limits. At each step, the maximum element in the range
is selected as a pivot, and visible connections are determined by scanning
outward while maintaining local maxima constraints.

This is an internal implementation used by `hvg`.
"""
function _hvg_core!(
    edges::Vector{Tuple{Int,Int}},
    x::AbstractVector{<:Real},
    left::Int,
    right::Int,
)
    # using a stack prevents StackOverflowError on massive time series
    stack = [(left, right)]

    while !isempty(stack)
        l, r = pop!(stack)
        r - l < 1 && continue

        # subarray maximum x[mid] is always visible to other nodes in the range
        mid = _argmax_range(x, l, r)

        # scan LEFT: connect mid to visible nodes in [l, mid-1]
        @inbounds if mid > l
            max_seen = x[mid - 1]
            push!(edges, (mid - 1, mid))
            for i in mid-2:-1:l
                xi = x[i]
                if xi > max_seen
                    push!(edges, (i, mid))
                    max_seen = xi
                end
            end
        end        

        # scan RIGHT: connect mid to visible nodes in [mid+1, r]
        @inbounds if mid < r
            max_seen = x[mid + 1]
            push!(edges, (mid, mid + 1))
            for j in mid+2:r
                xj = x[j]
                if xj > max_seen
                    push!(edges, (mid, j))
                    max_seen = xj
                end
            end
        end

        # recurse on independent sub-problems
        if l < mid - 1
            push!(stack, (l, mid - 1))
        end
        if mid + 1 < r
            push!(stack, (mid + 1, r))
        end
    end
end

"""
    _nvg_core!(edges, x, left, right)

Construct Natural Visibility Graph (NVG) edges for the subrange
`x[left:right]` and append them to `edges`.

The algorithm uses a divide-and-conquer strategy with an explicit stack
to avoid recursion limits. Visibility is determined using a geometric
criterion based on slope comparisons:

    (x[i] - x[mid]) / (mid - i)

is evaluated without division via cross-multiplication to improve numerical
stability and performance.

This is an internal implementation used by `nvg`.
"""
function _nvg_core!(
    edges::Vector{Tuple{Int,Int}},
    x::AbstractVector{<:Real},
    left::Int,
    right::Int,
)
    # statically infer the correct float type to prevent boxing/type-instability
    F = typeof(float(zero(eltype(x)))) 
    
    stack = [(left, right)]

    while !isempty(stack)
        l, r = pop!(stack)
        r - l < 1 && continue

        mid = _argmax_range(x, l, r)
        xm  = @inbounds x[mid]

        # scan LEFT
        # uses cross-multiplication: dy_new * dx_old > dy_old * dx_new
        # eliminates floating-point division from hot loop
        best_dy = typemin(F)
        best_dx = F(1) 

        @inbounds for i in mid-1:-1:l
            dy = F(x[i] - xm)
            dx = F(mid - i)
            if dy * best_dx > best_dy * dx
                push!(edges, (i, mid))
                best_dy = dy
                best_dx = dx
            end
        end

        # scan RIGHT
        best_dy = typemin(F)
        best_dx = F(1)

        @inbounds for j in mid+1:r
            dy = F(x[j] - xm)
            dx = F(j - mid)
            if dy * best_dx > best_dy * dx
                push!(edges, (mid, j))
                best_dy = dy
                best_dx = dx
            end
        end

        # recurse on independent sub-problems
        if l < mid - 1
            push!(stack, (l, mid - 1))
        end
        if mid + 1 < r
            push!(stack, (mid + 1, r))
        end
    end
end