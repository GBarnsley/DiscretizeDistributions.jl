using Test
using DiscretiseDistributions
using Distributions

@testset "DiscretiseDistributions.jl" begin
    include("test_discretise.jl")
    include("test_utility.jl")
end
