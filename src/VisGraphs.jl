module VisGraphs

using Plots: plot, plot!
using Plots.Colors: RGB
using Random: AbstractRNG, default_rng
using SparseArrays: sparse, issparse
using LinearAlgebra: Diagonal

include("utils/validation.jl")
include("core/internal.jl")
include("core/hvg.jl")
include("core/nvg.jl")
include("datasets/signals.jl")
include("analysis/metrics.jl")

export generate_sine, generate_random, generate_noisy_sine
export hvg, plot_hvg, nvg, plot_nvg
export whvg, plot_whvg, wnvg, plot_wnvg
export adjacency_matrix, degree_distribution, laplacian_matrix
export clustering_coefficient, average_path_length

end