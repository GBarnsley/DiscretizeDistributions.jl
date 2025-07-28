@testset "discretise function tests" begin
    
    @testset "Continuous distributions with fixed interval" begin
        # Test with Normal distribution
        normal_dist = Normal(0, 1)
        interval = 0.1
        discretised = discretise(normal_dist, interval)
        
        @test discretised isa DiscreteNonParametric
        @test sum(probs(discretised)) ≈ 1.0 atol=1e-10
        @test all(support(discretised) .>= minimum(support(discretised)))
        @test all(support(discretised) .<= maximum(support(discretised)))
        
        # Test that probabilities are positive
        @test all(probs(discretised) .> 0)
        
        # Test with custom quantiles
        discretised_custom = discretise(normal_dist, interval; min_quantile=0.01, max_quantile=0.99)
        @test discretised_custom isa DiscreteNonParametric
        @test sum(probs(discretised_custom)) ≈ 1.0 atol=1e-10
        
        # Test with Exponential distribution (unbounded above)
        exp_dist = Exponential(1.0)
        discretised_exp = discretise(exp_dist, 0.5)
        @test discretised_exp isa DiscreteNonParametric
        @test sum(probs(discretised_exp)) ≈ 1.0 atol=1e-10
        @test minimum(support(discretised_exp)) >= 0
        
        # Test with Uniform distribution (bounded)
        uniform_dist = Uniform(0, 10)
        discretised_uniform = discretise(uniform_dist, 1.0)
        @test discretised_uniform isa DiscreteNonParametric
        @test sum(probs(discretised_uniform)) ≈ 1.0 atol=1e-10
        @test minimum(support(discretised_uniform)) >= 0
        @test maximum(support(discretised_uniform)) <= 10
    end
    
    @testset "Continuous distributions with custom intervals" begin
        # Test with custom interval vector
        normal_dist = Normal(5, 2)
        custom_intervals = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0]
        discretised = discretise(normal_dist, custom_intervals)
        
        @test discretised isa DiscreteNonParametric
        @test sum(probs(discretised)) ≈ 1.0 atol=1e-10
        @test all(probs(discretised) .> 0)
        
        # Support should be based on the intervals (excluding last point)
        expected_support = custom_intervals[1:end-1]
        @test support(discretised) == expected_support[probs(discretised) .> 0]
        
        # Test with unsorted intervals (should be sorted internally)
        unsorted_intervals = [10.0, 0.0, 5.0, 2.0, 8.0]
        discretised_unsorted = discretise(normal_dist, unsorted_intervals)
        @test discretised_unsorted isa DiscreteNonParametric
        @test sum(probs(discretised_unsorted)) ≈ 1.0 atol=1e-10
    end
    
    @testset "Discrete distributions with fixed interval" begin
        # Test with Poisson distribution
        poisson_dist = Poisson(3.0)
        interval = 1
        discretised = discretise(poisson_dist, interval)
        
        @test discretised isa DiscreteNonParametric
        @test sum(probs(discretised)) ≈ 1.0 atol=1e-10
        @test all(probs(discretised) .> 0)
        @test minimum(support(discretised)) >= 0
        
        # Test with custom quantiles
        discretised_custom = discretise(poisson_dist, interval; min_quantile=0.01, max_quantile=0.95)
        @test discretised_custom isa DiscreteNonParametric
        @test sum(probs(discretised_custom)) ≈ 1.0 atol=1e-10
        
        # Test with Binomial distribution (bounded)
        binomial_dist = Binomial(10, 0.3)
        discretised_binomial = discretise(binomial_dist, 1)
        @test discretised_binomial isa DiscreteNonParametric
        @test sum(probs(discretised_binomial)) ≈ 1.0 atol=1e-10
        @test minimum(support(discretised_binomial)) >= 0
        @test maximum(support(discretised_binomial)) <= 10
        
        # Test with Geometric distribution (unbounded above)
        geometric_dist = Geometric(0.2)
        discretised_geom = discretise(geometric_dist, 2)
        @test discretised_geom isa DiscreteNonParametric
        @test sum(probs(discretised_geom)) ≈ 1.0 atol=1e-10
        @test minimum(support(discretised_geom)) >= 0
    end
    
    @testset "Discrete distributions with custom intervals" begin
        # Test with custom interval vector
        poisson_dist = Poisson(4.0)
        custom_intervals = [0, 2, 4, 6, 8, 12, 20]
        discretised = discretise(poisson_dist, custom_intervals)
        
        @test discretised isa DiscreteNonParametric
        @test sum(probs(discretised)) ≈ 1.0 atol=1e-10
        @test all(probs(discretised) .> 0)
        
        # Test with unsorted intervals
        unsorted_intervals = [8, 0, 4, 2, 12]
        discretised_unsorted = discretise(poisson_dist, unsorted_intervals)
        @test discretised_unsorted isa DiscreteNonParametric
        @test sum(probs(discretised_unsorted)) ≈ 1.0 atol=1e-10
    end
    
    @testset "Edge cases and error conditions" begin
        normal_dist = Normal(0, 1)
        
        # Test with very small interval
        small_interval = 0.001
        discretised_small = discretise(normal_dist, small_interval)
        @test discretised_small isa DiscreteNonParametric
        @test sum(probs(discretised_small)) ≈ 1.0 atol=1e-10
        
        # Test with very large interval
        large_interval = 10.0
        discretised_large = discretise(normal_dist, large_interval)
        @test discretised_large isa DiscreteNonParametric
        @test sum(probs(discretised_large)) ≈ 1.0 atol=1e-10
        
        # Test with single interval in vector
        single_interval = [0.0, 1.0]
        discretised_single = discretise(normal_dist, single_interval)
        @test discretised_single isa DiscreteNonParametric
        @test length(support(discretised_single)) <= 1
        
        # Test extreme quantiles
        discretised_extreme = discretise(normal_dist, 0.1; min_quantile=1e-6, max_quantile=1-1e-6)
        @test discretised_extreme isa DiscreteNonParametric
        @test sum(probs(discretised_extreme)) ≈ 1.0 atol=1e-10
    end
end
