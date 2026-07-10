using VisGraphs
using Test
using SparseArrays
using Random

@testset "VisGraphs.jl" begin
    include("signals_tests.jl")
    include("hvg_nvg_tests.jl")
    include("whvg_wnvg_tests.jl")
    include("analysis_tests.jl")
end