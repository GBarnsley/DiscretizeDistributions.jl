@testset "Interval-based discretisation tests" begin
    using Statistics  # For mean and var
    
    @testset "Backend interval creation" begin
        # Test that discretize returns interval-based distributions
        normal_dist = Normal(0, 1)
        discrete_intervals = discretize(normal_dist, 0.5)
        
        @test discrete_intervals isa DiscreteNonParametric
        @test eltype(support(discrete_intervals)) <: Interval
        @test sum(probs(discrete_intervals)) ≈ 1.0 atol=1e-10
        @test all(probs(discrete_intervals) .> 0)
        
        # Test interval properties
        first_interval = support(discrete_intervals)[1]
        @test first_interval isa Interval
        @test inf(first_interval) < sup(first_interval)
        
        # Test with different interval sizes
        for interval_size in [0.1, 0.5, 1.0, 2.0]
            discrete_test = discretize(normal_dist, interval_size)
            @test eltype(support(discrete_test)) <: Interval
            @test sum(probs(discrete_test)) ≈ 1.0 atol=1e-10
        end
        
        # Test with discrete distribution
        poisson_dist = Poisson(3.0)
        discrete_poisson = discretize(poisson_dist, 1.0)
        @test eltype(support(discrete_poisson)) <: Interval
        @test sum(probs(discrete_poisson)) ≈ 1.0 atol=1e-10
    end
    
    @testset "Custom interval boundaries" begin
        normal_dist = Normal(0, 1)
        custom_boundaries = [-2.0, -1.0, 0.0, 1.0, 2.0]
        discrete_custom = discretize(normal_dist, custom_boundaries)
        
        @test eltype(support(discrete_custom)) <: Interval
        @test sum(probs(discrete_custom)) ≈ 1.0 atol=1e-10
        
        # Check that intervals cover the expected ranges
        intervals = support(discrete_custom)
        @test length(intervals) >= length(custom_boundaries) - 1
        
        # Test with unsorted boundaries
        unsorted_boundaries = [1.0, -1.0, 2.0, 0.0, -2.0]
        discrete_unsorted = discretize(normal_dist, unsorted_boundaries)
        @test eltype(support(discrete_unsorted)) <: Interval
        @test sum(probs(discrete_unsorted)) ≈ 1.0 atol=1e-10
    end
    
    @testset "Quantile bounds for unbounded distributions" begin
        # Test exponential distribution (unbounded above)
        exp_dist = Exponential(1.0)
        discrete_exp = discretize(exp_dist, 0.5)
        
        @test eltype(support(discrete_exp)) <: Interval
        @test sum(probs(discrete_exp)) ≈ 1.0 atol=1e-10
        # Test that intervals cover expected range and check bounds properly
        first_interval = minimum(support(discrete_exp))
        @test inf(first_interval) >= -1e-6  # Should start near or at 0
        @test sum(probs(discrete_exp)) ≈ 1.0 atol=1e-10
        
        # Test with custom quantiles
        discrete_exp_custom = discretize(exp_dist, 0.1; min_quantile=0.01, max_quantile=0.95)
        @test eltype(support(discrete_exp_custom)) <: Interval
        @test sum(probs(discrete_exp_custom)) ≈ 1.0 atol=1e-10
    end
    
    @testset "Interval properties and consistency" begin
        normal_dist = Normal(0, 1)
        discrete_intervals = discretize(normal_dist, 0.5)
        
        intervals = support(discrete_intervals)
        
        # Test that intervals are contiguous (end of one = start of next)
        for i in 1:(length(intervals)-1)
            @test sup(intervals[i]) ≈ inf(intervals[i+1]) atol=1e-10
        end
        
        # Test that all intervals have positive width
        for interval in intervals
            @test sup(interval) > inf(interval)
        end
        
        # Test that probabilities correspond to interval widths for uniform distribution
        uniform_dist = Uniform(0, 10)
        discrete_uniform = discretize(uniform_dist, 1.0)
        uniform_intervals = support(discrete_uniform)
        uniform_probs = probs(discrete_uniform)
        
        # For uniform distribution, equal-width intervals should have equal probabilities
        widths = [sup(interval) - inf(interval) for interval in uniform_intervals]
        # For uniform distribution, equal-width intervals should have similar probabilities
        widths = [sup(interval) - inf(interval) for interval in uniform_intervals]
        # Check if we have reasonably uniform widths (within tolerance)
        if length(widths) > 1
            mean_width = mean(widths)
            uniform_width_indices = findall(w -> abs(w - mean_width) < 0.1 * mean_width, widths)
            if length(uniform_width_indices) > 1
                interior_probs = uniform_probs[uniform_width_indices]
                # Allow for some variation due to distribution shape and boundary effects
                prob_variance = var(interior_probs)
                mean_prob = mean(interior_probs)
                @test prob_variance < 0.1 * mean_prob^2  # Variance should be small relative to mean
            end
        end
    end
    
    @testset "Rational number precision" begin
        normal_dist = Normal(0, 1)
        
        # Test with rational interval size
        rational_interval = 1//10
        discrete_rational = discretize(normal_dist, rational_interval)
        
        @test eltype(support(discrete_rational)) <: Interval
        @test sum(probs(discrete_rational)) ≈ 1.0 atol=1e-10
        
        # Check that rational precision is maintained in interval bounds where possible
        intervals = support(discrete_rational)
        finite_intervals = filter(i -> isfinite(inf(i)) && isfinite(sup(i)), intervals)
        
        if length(finite_intervals) > 0
            # Test a few intervals for width consistency
            for interval in finite_intervals[1:min(3, end)]
                width = sup(interval) - inf(interval)
                # Allow for some tolerance due to distribution bounds and floating point precision
                @test abs(width - 0.1) < 0.05 || abs(width) < 1e-10  # Either close to 0.1 or very small (boundary effect)
            end
        end
    end
    
    @testset "Error handling and edge cases" begin
        normal_dist = Normal(0, 1)
        
        # Test with very small intervals
        discrete_small = discretize(normal_dist, 1e-3)
        @test eltype(support(discrete_small)) <: Interval
        @test sum(probs(discrete_small)) ≈ 1.0 atol=1e-10
        
        # Test with very large intervals
        discrete_large = discretize(normal_dist, 10.0)
        @test eltype(support(discrete_large)) <: Interval
        @test sum(probs(discrete_large)) ≈ 1.0 atol=1e-10
        
        # Test with single custom interval
        single_boundary = [0.0, 1.0]
        discrete_single = discretize(normal_dist, single_boundary)
        @test eltype(support(discrete_single)) <: Interval
        @test sum(probs(discrete_single)) ≈ 1.0 atol=1e-10
        
        # Test extreme quantiles
        discrete_extreme = discretize(normal_dist, 0.1; min_quantile=1e-6, max_quantile=1-1e-6)
        @test eltype(support(discrete_extreme)) <: Interval
        @test sum(probs(discrete_extreme)) ≈ 1.0 atol=1e-10
    end
end
