using Test
using DiscretiseDistributions
using Distributions, IntervalArithmetic
@testset "DiscretiseDistributions.jl" begin
    include("test_discretise.jl")
    include("test_utility.jl")
    include("test_intervals.jl")
end
