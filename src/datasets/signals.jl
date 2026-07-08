"""
    generate_sine(n::Integer=100)

Generate a sine wave sampled at `n` evenly spaced points over the interval
``[0, 4π]``.

The returned signal spans two complete periods of the sine wave (each
period is `2π`, and the sampling interval is `4π`), and is primarily
intended for testing, examples, and benchmarking algorithms that operate
on time series.

# Examples
```jldoctest
julia> x = generate_sine(5);

julia> length(x)
5

julia> x[1] ≈ 0.0
true
```
"""
function generate_sine(n::Integer=100)
    return sin.(range(0, 4π, length=n))
end

"""
    generate_random(n::Integer=100)

Generate a random time series of length `n`.

Each sample is drawn independently from a uniform distribution on the interval
``[0, 1)`` using Julia's default random number generator. The resulting series
is useful for testing and benchmarking algorithms on unstructured data.

# Examples
```jldoctest
julia> x = generate_random(10);

julia> length(x)
10

julia> all(0 .<= generate_random(10) .< 1)
true
```
"""
function generate_random(n::Integer=100)
    return rand(n)
end


"""
    generate_noisy_sine(n::Integer=100, noise::Real=0.2)

Generate a sine wave with additive Gaussian noise.

The underlying signal is sampled at `n` evenly spaced points over the interval
``[0, 4π]``. Independent Gaussian noise with standard deviation `noise` is
added to each sample.

This function is useful for testing the robustness of time-series algorithms
under noisy conditions.

# Examples
```jldoctest
julia> x = generate_noisy_sine(3, 0.0);

julia> y = generate_sine(3);

julia> x == y
true
```
"""
function generate_noisy_sine(n::Integer=100, noise::Real=0.2)
    return sin.(range(0, 4π, length=n)) .+ noise .* randn(n)
end
