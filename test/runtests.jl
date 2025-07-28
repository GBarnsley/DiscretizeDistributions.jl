using Test
using DiscretizeDistributions
using Distributions, IntervalArithmetic
@testset "DiscretizeDistributions.jl" begin
    include("test_discretize.jl")
    include("test_utility.jl")
    include("test_intervals.jl")
end
