"""
    _validate(x)

Check that `x` is a valid input time series.

Ensures that:
- `x` contains at least two elements
- all elements are finite real numbers (no `NaN` or `Inf`)
- `x` has standard 1-based linear indexing (no `OffsetArray` or custom indices)

Throws an `ArgumentError` if validation fails.

This is an internal precondition check used by all visibility graph
construction algorithms.
"""
function _validate(x::AbstractVector{<:Real})
    Base.require_one_based_indexing(x)
    length(x) ≥ 2 || throw(ArgumentError("`x` must contain at least two values."))
    all(isfinite, x) || throw(ArgumentError("`x` must contain only finite values."))
end

"""
    _validate_length(n)

Check that `n` is a valid time series length.

Ensures that:
- `n` is at least two

Throws an `ArgumentError` if validation fails.

This is an internal precondition check used by signal generation
functions that require a minimum number of samples.
"""
function _validate_length(n::Integer)
    n ≥ 2 || throw(ArgumentError("`n` must be at least two."))
end