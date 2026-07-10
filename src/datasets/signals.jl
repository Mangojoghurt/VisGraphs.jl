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
    _validate_length(n)
    return sin.(range(0, 4π, length=n))
end

"""
    generate_random(n::Integer=100; rng::AbstractRNG=default_rng())

Generate a random time series of length `n`.

Each sample is drawn independently from a uniform distribution on the interval
``[0, 1)`` using the provided random number generator `rng`. By default,
Julia's default random number generator is used.

The `rng` argument can be used to obtain reproducible results by passing a
seeded random number generator. The resulting series is useful for testing and
benchmarking algorithms on unstructured data.

# Examples
```jldoctest
julia> x = generate_random(10);

julia> length(x)
10

julia> all(0 .<= generate_random(10) .< 1)
true
```
"""
function generate_random(n::Integer=100; rng::AbstractRNG=default_rng())
    _validate_length(n)
    return rand(rng, n)
end


"""
    generate_noisy_sine(n::Integer=100, noise::Real=0.2; rng::AbstractRNG=default_rng())

Generate a sine wave with additive Gaussian noise.

The underlying signal is sampled at `n` evenly spaced points over the interval
``[0, 4π]``. Independent Gaussian noise with standard deviation `noise` is
added to each sample using the provided random number generator `rng`.

By default, Julia's default random number generator is used. The `rng`
argument can be used to obtain reproducible noisy signals by passing a seeded
random number generator.

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
function generate_noisy_sine(n::Integer=100, noise::Real=0.2; rng::AbstractRNG=default_rng())
    _validate_length(n)
    return sin.(range(0, 4π, length=n)) .+ noise .* randn(rng, n)
end
