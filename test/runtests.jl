using Test
using discretizeDistributions
using Distributions, IntervalArithmetic
@testset "discretizeDistributions.jl" begin
    include("test_discretize.jl")
    include("test_utility.jl")
    include("test_intervals.jl")
end
