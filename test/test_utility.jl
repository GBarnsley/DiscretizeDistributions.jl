@testset "Utility function tests" begin
    
    @testset "center_distribution function tests" begin
        # Create a simple discrete distribution for testing
        support_vals = [1.0, 2.0, 3.0, 4.0, 5.0]
        probabilities = [0.1, 0.2, 0.4, 0.2, 0.1]
        test_dist = DiscreteNonParametric(support_vals, probabilities)
        
        # Test center_distribution without interval parameter
        centered = center_distribution(test_dist)
        @test centered isa DiscreteNonParametric
        
        # Check that support is centered between original intervals
        original_support = support(test_dist)
        centered_support = support(centered)
        expected_centers = original_support[1:end-1] .+ (diff(original_support) ./ 2)
        @test centered_support == expected_centers
        
        # Check that probabilities are preserved (except last one is dropped)
        @test probs(centered) == probs(test_dist)[1:end-1]
        @test sum(probs(centered)) ≈ sum(probabilities[1:end-1]) atol=1e-10
        
        # Test with uniform spacing
        uniform_support = [0.0, 1.0, 2.0, 3.0, 4.0]
        uniform_probs = [0.25, 0.25, 0.25, 0.25, 0.0]  # Last probability is 0
        uniform_dist = DiscreteNonParametric(uniform_support, uniform_probs)
        centered_uniform = center_distribution(uniform_dist)
        
        expected_uniform_centers = [0.5, 1.5, 2.5, 3.5]
        @test support(centered_uniform) == expected_uniform_centers
        
        # Test center_distribution with interval parameter
        interval = 0.5
        centered_with_interval = center_distribution(test_dist, interval)
        @test centered_with_interval isa DiscreteNonParametric
        
        # Check that support is shifted by interval/2
        expected_shifted = support(test_dist) .+ (interval / 2)
        @test support(centered_with_interval) == expected_shifted
        
        # Check that probabilities are preserved
        @test probs(centered_with_interval) == probs(test_dist)
        @test sum(probs(centered_with_interval)) ≈ 1.0 atol=1e-10
        
        # Test with different interval values
        for test_interval in [0.1, 1.0, 2.5, 10.0]
            centered_test = center_distribution(test_dist, test_interval)
            expected_support = support(test_dist) .+ (test_interval / 2)
            @test support(centered_test) == expected_support
            @test probs(centered_test) == probs(test_dist)
        end
        
        # Test with negative interval
        negative_interval = -0.3
        centered_negative = center_distribution(test_dist, negative_interval)
        expected_negative = support(test_dist) .+ (negative_interval / 2)
        @test support(centered_negative) == expected_negative
    end
    
    @testset "right_align_distribution function tests" begin
        # Create a simple discrete distribution for testing
        support_vals = [1.0, 2.0, 3.0, 4.0, 5.0]
        probabilities = [0.1, 0.2, 0.4, 0.2, 0.1]
        test_dist = DiscreteNonParametric(support_vals, probabilities)
        
        # Test right_align_distribution without interval parameter
        right_aligned = right_align_distribution(test_dist)
        @test right_aligned isa DiscreteNonParametric
        
        # Check that support is shifted to the right (using second through last elements)
        original_support = support(test_dist)
        aligned_support = support(right_aligned)
        expected_right = original_support[2:end]
        @test aligned_support == expected_right
        
        # Check that probabilities are from first to second-to-last
        @test probs(right_aligned) == probs(test_dist)[1:end-1]
        @test sum(probs(right_aligned)) ≈ sum(probabilities[1:end-1]) atol=1e-10
        
        # Test with uniform spacing
        uniform_support = [0.0, 1.0, 2.0, 3.0, 4.0]
        uniform_probs = [0.25, 0.25, 0.25, 0.25, 0.0]  # Last probability is 0
        uniform_dist = DiscreteNonParametric(uniform_support, uniform_probs)
        right_aligned_uniform = right_align_distribution(uniform_dist)
        
        expected_uniform_right = [1.0, 2.0, 3.0, 4.0]
        @test support(right_aligned_uniform) == expected_uniform_right
        
        # Test right_align_distribution with interval parameter
        interval = 0.7
        right_aligned_with_interval = right_align_distribution(test_dist, interval)
        @test right_aligned_with_interval isa DiscreteNonParametric
        
        # Check that support is shifted by interval
        expected_shifted = support(test_dist) .+ interval
        @test support(right_aligned_with_interval) == expected_shifted
        
        # Check that probabilities are preserved
        @test probs(right_aligned_with_interval) == probs(test_dist)
        @test sum(probs(right_aligned_with_interval)) ≈ 1.0 atol=1e-10
        
        # Test with different interval values
        for test_interval in [0.1, 1.0, 2.5, 10.0, -1.5]
            aligned_test = right_align_distribution(test_dist, test_interval)
            expected_support = support(test_dist) .+ test_interval
            @test support(aligned_test) == expected_support
            @test probs(aligned_test) == probs(test_dist)
        end
        
        # Test with zero interval
        zero_interval = 0.0
        aligned_zero = right_align_distribution(test_dist, zero_interval)
        @test support(aligned_zero) == support(test_dist)
        @test probs(aligned_zero) == probs(test_dist)
    end
    
    @testset "Edge cases for utility functions" begin
        # Test with single-point distribution
        single_support = [5.0]
        single_probs = [1.0]
        single_dist = DiscreteNonParametric(single_support, single_probs)
        
        # center_distribution with single point should error
        @test_throws ErrorException center_distribution(single_dist)
        
        # right_align_distribution with single point should error
        @test_throws ErrorException right_align_distribution(single_dist)

        # But with interval parameter, should preserve the distribution
        centered_single_interval = center_distribution(single_dist, 1.0)
        @test support(centered_single_interval) == [5.5]
        @test probs(centered_single_interval) == [1.0]
        
        right_aligned_single_interval = right_align_distribution(single_dist, 1.0)
        @test support(right_aligned_single_interval) == [6.0]
        @test probs(right_aligned_single_interval) == [1.0]
        
        # Test with two-point distribution
        two_support = [1.0, 3.0]
        two_probs = [0.6, 0.4]
        two_dist = DiscreteNonParametric(two_support, two_probs)
        
        centered_two = center_distribution(two_dist)
        @test support(centered_two) == [2.0]  # (1 + 3)/2 = 2, but only first prob
        @test probs(centered_two) == [0.6]
        
        right_aligned_two = right_align_distribution(two_dist)
        @test support(right_aligned_two) == [3.0]  # second support value
        @test probs(right_aligned_two) == [0.6]    # first probability
    end
end
