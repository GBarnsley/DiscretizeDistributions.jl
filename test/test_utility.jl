@testset "Utility function tests for interval-based backend" begin
    
    # Import the internal function for testing
    import DiscretiseDistributions: remove_infinities
    
    @testset "remove_infinities function tests" begin
        # Create a distribution using discretise function then manually add infinite intervals for testing
        normal_dist = Normal(0, 1)
        base_discrete = discretise(normal_dist, 0.5)
        
        # We'll test remove_infinities by creating a distribution that already includes infinite intervals
        # through the discretise function with extreme quantiles
        discrete_with_extremes = discretise(normal_dist, 0.5; min_quantile=1e-10, max_quantile=1-1e-10)
        
        # Test remove_infinities on this distribution
        cleaned = remove_infinities(discrete_with_extremes)
        @test cleaned isa DiscreteNonParametric
        @test eltype(support(cleaned)) <: Interval
        
        # Check that infinite intervals are removed (if any existed)
        cleaned_intervals = support(cleaned)
        @test all(interval -> !isinf(inf(interval)) && !isinf(sup(interval)), cleaned_intervals)
        
        # Check that probabilities are renormalized
        @test sum(probs(cleaned)) ≈ 1.0 atol=1e-10
        
        # Test with a distribution that has finite intervals (should not change)
        finite_dist = discretise(Normal(0, 1), 0.5; min_quantile=0.01, max_quantile=0.99)
        cleaned_finite = remove_infinities(finite_dist)
        @test length(support(cleaned_finite)) <= length(support(finite_dist))
        @test sum(probs(cleaned_finite)) ≈ 1.0 atol=1e-10
    end
    
    @testset "left_align_distribution function tests" begin
        # Create interval-based distribution using discretise function
        normal_dist = Normal(0, 1)
        interval_dist = discretise(normal_dist, 0.5; min_quantile=0.01, max_quantile=0.99)
        
        # Test left alignment
        left_aligned = left_align_distribution(interval_dist)
        @test left_aligned isa DiscreteNonParametric
        @test eltype(support(left_aligned)) <: Real
        
        # Check that support uses left bounds (inf) of intervals (after removing infinities)
        # The alignment functions correctly remove infinite intervals
        finite_intervals = [interval for interval in support(interval_dist) if isfinite(inf(interval)) && isfinite(sup(interval))]
        expected_support = [inf(interval) for interval in finite_intervals]
        @test support(left_aligned) == expected_support
        @test length(probs(left_aligned)) == length(finite_intervals)  # Should match finite intervals only
        @test sum(probs(left_aligned)) ≈ 1.0 atol=1e-10
        
        # Test with distribution that has been cleaned of infinities
        infinite_dist = discretise(normal_dist, 0.5; min_quantile=1e-10, max_quantile=1-1e-10)
        cleaned_and_aligned = left_align_distribution(remove_infinities(infinite_dist))
        @test eltype(support(cleaned_and_aligned)) <: Real
        @test all(isfinite, support(cleaned_and_aligned))
        @test sum(probs(cleaned_and_aligned)) ≈ 1.0 atol=1e-10
        
        # Test with a more controlled example
        uniform_dist = Uniform(0, 10)
        uniform_intervals = discretise(uniform_dist, 1.0)
        left_uniform = left_align_distribution(uniform_intervals)
        @test eltype(support(left_uniform)) <: Real
        @test all(isfinite, support(left_uniform))
    end
    
    @testset "centred_distribution function tests" begin
        # Create interval-based distribution using discretise function
        normal_dist = Normal(0, 1)
        interval_dist = discretise(normal_dist, 0.5; min_quantile=0.01, max_quantile=0.99)
        
        # Test centred alignment
        centred = centred_distribution(interval_dist)
        @test centred isa DiscreteNonParametric
        @test eltype(support(centred)) <: Real
        
        # Check that support uses midpoints of intervals (after removing infinities)
        finite_intervals = [interval for interval in support(interval_dist) if isfinite(inf(interval)) && isfinite(sup(interval))]
        expected_support = [(inf(interval) + sup(interval)) / 2 for interval in finite_intervals]
        @test support(centred) ≈ expected_support atol=1e-10
        @test length(probs(centred)) == length(finite_intervals)  # Should match finite intervals only
        @test sum(probs(centred)) ≈ 1.0 atol=1e-10
        
        # Test with uniform intervals for more predictable midpoints
        uniform_dist = Uniform(0, 10)
        uniform_intervals = discretise(uniform_dist, 1.0)
        centred_uniform = centred_distribution(uniform_intervals)
        
        # Check that midpoints are correctly calculated
        uniform_test_intervals = support(uniform_intervals)
        expected_uniform = [(inf(interval) + sup(interval)) / 2 for interval in uniform_test_intervals]
        @test support(centred_uniform) ≈ expected_uniform atol=1e-10
        
        # Test edge case with narrow distribution
        narrow_dist = Normal(5, 0.1)  # Very narrow normal distribution
        narrow_intervals = discretise(narrow_dist, 0.1)
        centred_narrow = centred_distribution(narrow_intervals)
        @test eltype(support(centred_narrow)) <: Real
        @test all(isfinite, support(centred_narrow))
    end
    
    @testset "right_align_distribution function tests" begin
        # Create interval-based distribution using discretise function  
        normal_dist = Normal(0, 1)
        interval_dist = discretise(normal_dist, 0.5; min_quantile=0.01, max_quantile=0.99)
        
        # Test right alignment
        right_aligned = right_align_distribution(interval_dist)
        @test right_aligned isa DiscreteNonParametric
        @test eltype(support(right_aligned)) <: Real
        
        # Check that support uses right bounds (sup) of intervals (after removing infinities)
        finite_intervals = [interval for interval in support(interval_dist) if isfinite(inf(interval)) && isfinite(sup(interval))]
        expected_support = [sup(interval) for interval in finite_intervals]
        @test support(right_aligned) == expected_support
        @test length(probs(right_aligned)) == length(finite_intervals)  # Should match finite intervals only
        @test sum(probs(right_aligned)) ≈ 1.0 atol=1e-10
        
        # Test with uniform distribution for more predictable behavior
        uniform_dist = Uniform(0, 10)
        uniform_intervals = discretise(uniform_dist, 1.0)
        right_uniform = right_align_distribution(uniform_intervals)
        
        # Check that all support values are finite and properly aligned
        @test eltype(support(right_uniform)) <: Real
        @test all(isfinite, support(right_uniform))
        @test sum(probs(right_uniform)) ≈ 1.0 atol=1e-10
        
        # Verify that right bounds are used
        uniform_test_intervals = support(uniform_intervals)
        expected_uniform = [sup(interval) for interval in uniform_test_intervals]
        @test support(right_uniform) == expected_uniform
    end
    
    @testset "Integration tests with discretise function" begin
        # Test the complete workflow: discretise -> align
        normal_dist = Normal(0, 1)
        discretised = discretise(normal_dist, 0.5; min_quantile=0.01, max_quantile=0.99)
        
        # Test that we can apply all alignment functions
        left_aligned = left_align_distribution(discretised)
        centred = centred_distribution(discretised)
        right_aligned = right_align_distribution(discretised)
        
        @test eltype(support(left_aligned)) <: Real
        @test eltype(support(centred)) <: Real  
        @test eltype(support(right_aligned)) <: Real
        
        @test sum(probs(left_aligned)) ≈ 1.0 atol=1e-10
        @test sum(probs(centred)) ≈ 1.0 atol=1e-10
        @test sum(probs(right_aligned)) ≈ 1.0 atol=1e-10
        
        # Test that alignment preserves probability mass for cleaned distributions
        cleaned = remove_infinities(discretised)
        cleaned_left = left_align_distribution(cleaned)
        cleaned_centred = centred_distribution(cleaned)
        cleaned_right = right_align_distribution(cleaned)
        
        @test length(probs(cleaned_left)) == length(probs(cleaned))
        @test length(probs(cleaned_centred)) == length(probs(cleaned))
        @test length(probs(cleaned_right)) == length(probs(cleaned))
        
        # Test relationship between alignments for uniform intervals
        uniform_dist = Uniform(0, 10)
        uniform_discretised = discretise(uniform_dist, 1.0)
        uniform_intervals = support(uniform_discretised)
        
        if length(uniform_intervals) > 1
            uniform_left = left_align_distribution(uniform_discretised)
            uniform_right = right_align_distribution(uniform_discretised)
            uniform_center = centred_distribution(uniform_discretised)
            
            # Check that intervals have consistent structure for finite intervals
            left_vals = support(uniform_left)
            right_vals = support(uniform_right)
            center_vals = support(uniform_center)
            
            # Test with finite intervals only
            finite_indices = findall(i -> isfinite(inf(uniform_intervals[i])) && isfinite(sup(uniform_intervals[i])), 1:length(uniform_intervals))
            
            for i in finite_indices[1:min(3, end)]  # Test first few finite intervals
                interval_width = sup(uniform_intervals[i]) - inf(uniform_intervals[i])
                if isfinite(interval_width) && interval_width > 0
                    @test right_vals[i] ≈ left_vals[i] + interval_width atol=1e-10
                    @test center_vals[i] ≈ left_vals[i] + interval_width/2 atol=1e-10
                end
            end
        end
        
        # Test with removal of infinities
        cleaned = remove_infinities(discretised)
        cleaned_left = left_align_distribution(cleaned)
        cleaned_centred = centred_distribution(cleaned)
        cleaned_right = right_align_distribution(cleaned)
        
        @test all(isfinite, support(cleaned_left))
        @test all(isfinite, support(cleaned_centred))
        @test all(isfinite, support(cleaned_right))
    end
    
    @testset "Edge cases and error conditions" begin
        # Test with single interval
        single_dist = discretise(Normal(0, 0.1), 2.0)  # Wide interval on narrow distribution
        
        if length(support(single_dist)) >= 1
            left_single = left_align_distribution(single_dist)
            centred_single = centred_distribution(single_dist)
            right_single = right_align_distribution(single_dist)
            
            @test eltype(support(left_single)) <: Real
            @test eltype(support(centred_single)) <: Real
            @test eltype(support(right_single)) <: Real
            
            @test sum(probs(left_single)) ≈ 1.0 atol=1e-10
            @test sum(probs(centred_single)) ≈ 1.0 atol=1e-10
            @test sum(probs(right_single)) ≈ 1.0 atol=1e-10
        end
        
        # Test with distribution that creates extreme values
        extreme_dist = discretise(Normal(0, 1), 0.5; min_quantile=1e-10, max_quantile=1-1e-10)
        extreme_cleaned = remove_infinities(extreme_dist)
        
        if length(support(extreme_cleaned)) > 0
            @test all(isfinite, support(left_align_distribution(extreme_cleaned)))
            @test all(isfinite, support(centred_distribution(extreme_cleaned)))
            @test all(isfinite, support(right_align_distribution(extreme_cleaned)))
        end
        
        # Test with very small intervals to check numerical stability
        tiny_dist = discretise(Normal(0, 1), 0.001; min_quantile=0.01, max_quantile=0.99)
        tiny_cleaned = remove_infinities(tiny_dist)
        
        tiny_centred = centred_distribution(tiny_cleaned)
        @test all(isfinite, support(tiny_centred))
        @test sum(probs(tiny_centred)) ≈ 1.0 atol=1e-10
    end
end
