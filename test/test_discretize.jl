@testset "discretize function tests with interval backend" begin
    @testset "Continuous distributions with fixed interval" begin
        # Test with Normal distribution
        normal_dist = Normal(0, 1)
        interval = 0.1
        discretized = discretize(normal_dist, interval)

        @test discretized isa DiscreteNonParametric
        @test eltype(support(discretized)) <: Interval
        @test sum(probs(discretized))≈1.0 atol=1e-10

        # Test that probabilities are positive
        @test all(probs(discretized) .> 0)

        # Test interval properties
        intervals = support(discretized)
        @test all(interval -> inf(interval) < sup(interval), intervals)

        # Test with custom quantiles
        discretized_custom = discretize(
            normal_dist, interval; min_quantile = 0.01, max_quantile = 0.99)
        @test discretized_custom isa DiscreteNonParametric
        @test eltype(support(discretized_custom)) <: Interval
        @test sum(probs(discretized_custom))≈1.0 atol=1e-10

        # Test with Exponential distribution (unbounded above)
        exp_dist = Exponential(1.0)
        discretized_exp = discretize(exp_dist, 0.5)
        @test discretized_exp isa DiscreteNonParametric
        @test eltype(support(discretized_exp)) <: Interval
        @test sum(probs(discretized_exp))≈1.0 atol=1e-10

        # Check first interval starts at or near 0
        first_interval = minimum(support(discretized_exp))
        @test inf(first_interval) >= -1e-10  # Allow for small numerical errors

        # Test with Uniform distribution (bounded)
        uniform_dist = Uniform(0, 10)
        discretized_uniform = discretize(uniform_dist, 1.0)
        @test discretized_uniform isa DiscreteNonParametric
        @test eltype(support(discretized_uniform)) <: Interval
        @test sum(probs(discretized_uniform))≈1.0 atol=1e-10

        # Check interval bounds cover the distribution support
        first_interval = minimum(support(discretized_uniform))
        last_interval = maximum(support(discretized_uniform))
        @test inf(first_interval) >= -1e-10
        @test sup(last_interval) <= 10.0 + 1e-10
    end

    @testset "Continuous distributions with custom intervals" begin
        # Test with custom interval vector
        normal_dist = Normal(5, 2)
        custom_intervals = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0]
        discretized = discretize(normal_dist, custom_intervals)

        @test discretized isa DiscreteNonParametric
        @test eltype(support(discretized)) <: Interval
        @test sum(probs(discretized))≈1.0 atol=1e-10
        @test all(probs(discretized) .> 0)

        # Test interval structure corresponds to boundaries
        intervals = support(discretized)
        @test length(intervals) >= 1  # Should have at least one interval
        @test length(intervals) <= length(custom_intervals) + 2  # May include infinite tail intervals

        # Test with unsorted intervals (should be sorted internally)
        unsorted_intervals = [10.0, 0.0, 5.0, 2.0, 8.0]
        discretized_unsorted = discretize(normal_dist, unsorted_intervals)
        @test discretized_unsorted isa DiscreteNonParametric
        @test eltype(support(discretized_unsorted)) <: Interval
        @test sum(probs(discretized_unsorted))≈1.0 atol=1e-10
    end

    @testset "Discrete distributions with fixed interval" begin
        # Test with Poisson distribution
        poisson_dist = Poisson(3.0)
        interval = 1
        discretized = discretize(poisson_dist, interval)

        @test discretized isa DiscreteNonParametric
        @test eltype(support(discretized)) <: Interval
        @test sum(probs(discretized))≈1.0 atol=1e-10
        @test all(probs(discretized) .> 0)

        # Check first interval starts at or near 0
        first_interval = minimum(support(discretized))
        @test inf(first_interval) >= -1e-10

        # Test with custom quantiles
        discretized_custom = discretize(
            poisson_dist, interval; min_quantile = 0.01, max_quantile = 0.95)
        @test discretized_custom isa DiscreteNonParametric
        @test eltype(support(discretized_custom)) <: Interval
        @test sum(probs(discretized_custom))≈1.0 atol=1e-10

        # Test with Binomial distribution (bounded)
        binomial_dist = Binomial(10, 0.3)
        discretized_binomial = discretize(binomial_dist, 1)
        @test discretized_binomial isa DiscreteNonParametric
        @test eltype(support(discretized_binomial)) <: Interval
        @test sum(probs(discretized_binomial))≈1.0 atol=1e-10

        # Check interval bounds are reasonable for bounded distribution
        first_interval = minimum(support(discretized_binomial))
        last_interval = maximum(support(discretized_binomial))
        @test inf(first_interval) >= -1e-10
        @test sup(last_interval) <= 10.0 + 1e-10

        # Test with Geometric distribution (unbounded above)
        geometric_dist = Geometric(0.2)
        discretized_geom = discretize(geometric_dist, 2)
        @test discretized_geom isa DiscreteNonParametric
        @test eltype(support(discretized_geom)) <: Interval
        @test sum(probs(discretized_geom))≈1.0 atol=1e-10

        # Check first interval starts at or near 0
        first_interval = minimum(support(discretized_geom))
        @test inf(first_interval) >= -1e-10
    end

    @testset "Discrete distributions with custom intervals" begin
        # Test with custom interval vector
        poisson_dist = Poisson(4.0)
        custom_intervals = [0, 2, 4, 6, 8, 12, 20]
        discretized = discretize(poisson_dist, custom_intervals)

        @test discretized isa DiscreteNonParametric
        @test eltype(support(discretized)) <: Interval
        @test sum(probs(discretized))≈1.0 atol=1e-10
        @test all(probs(discretized) .> 0)

        # Test with unsorted intervals
        unsorted_intervals = [8, 0, 4, 2, 12]
        discretized_unsorted = discretize(poisson_dist, unsorted_intervals)
        @test discretized_unsorted isa DiscreteNonParametric
        @test eltype(support(discretized_unsorted)) <: Interval
        @test sum(probs(discretized_unsorted))≈1.0 atol=1e-10
    end

    @testset "Edge cases and error conditions" begin
        normal_dist = Normal(0, 1)

        # Test with very small interval
        small_interval = 0.001
        discretized_small = discretize(normal_dist, small_interval)
        @test discretized_small isa DiscreteNonParametric
        @test eltype(support(discretized_small)) <: Interval
        @test sum(probs(discretized_small))≈1.0 atol=1e-10

        # Test with very large interval
        large_interval = 10.0
        discretized_large = discretize(normal_dist, large_interval)
        @test discretized_large isa DiscreteNonParametric
        @test eltype(support(discretized_large)) <: Interval
        @test sum(probs(discretized_large))≈1.0 atol=1e-10

        # Test with single interval in vector
        single_interval = [0.0, 1.0]
        discretized_single = discretize(normal_dist, single_interval)
        @test discretized_single isa DiscreteNonParametric
        @test eltype(support(discretized_single)) <: Interval
        @test length(support(discretized_single)) >= 1  # Should have at least one interval

        # Test extreme quantiles
        discretized_extreme = discretize(
            normal_dist, 0.1; min_quantile = 1e-6, max_quantile = 1 - 1e-6)
        @test discretized_extreme isa DiscreteNonParametric
        @test eltype(support(discretized_extreme)) <: Interval
        @test sum(probs(discretized_extreme))≈1.0 atol=1e-10
    end
end
