using Test
using DiscretizeDistributions
using Distributions, IntervalArithmetic
using Aqua

# Check if JET is available
jet_available = try
    using JET
    true
catch
    false
end

@testset "DiscretizeDistributions.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(DiscretizeDistributions)
    end

    if jet_available
        @testset "JET" begin
            JET.test_package(DiscretizeDistributions)
        end
    end

    include("test_discretize.jl")
    include("test_utility.jl")
    include("test_intervals.jl")
end
