using Test
using DiscretizeDistributions
using Distributions, IntervalArithmetic
using Aqua, JET

@testset "DiscretizeDistributions.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(DiscretizeDistributions)
    end

    @testset "JET" begin
        JET.test_package(DiscretizeDistributions)
    end

    include("test_discretize.jl")
    include("test_utility.jl")
    include("test_intervals.jl")
end
