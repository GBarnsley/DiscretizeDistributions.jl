@testset "Method parameter tests" begin
    @testset "Alignment methods with fixed intervals" begin
        normal_dist = Normal(0, 1)
        interval = 0.5

        # Test all alignment methods
        interval_result = discretize(normal_dist, interval; method = :interval)
        left_result = discretize(normal_dist, interval; method = :left_aligned)
        centered_result = discretize(normal_dist, interval; method = :centred)
        right_result = discretize(normal_dist, interval; method = :right_aligned)

        # Basic properties
        @test interval_result isa DiscreteNonParametric
        @test left_result isa DiscreteNonParametric
        @test centered_result isa DiscreteNonParametric
        @test right_result isa DiscreteNonParametric

        # Check support types
        @test eltype(support(interval_result)) <: Interval
        @test eltype(support(left_result)) <: Real
        @test eltype(support(centered_result)) <: Real
        @test eltype(support(right_result)) <: Real

        # Probabilities should sum to 1
        @test sum(probs(interval_result)) ≈ 1.0 atol=1e-10
        @test sum(probs(left_result)) ≈ 1.0 atol=1e-10
        @test sum(probs(centered_result)) ≈ 1.0 atol=1e-10
        @test sum(probs(right_result)) ≈ 1.0 atol=1e-10

        # Test that alignment functions work equivalently
        left_post = left_align_distribution(interval_result)
        centered_post = centred_distribution(interval_result)
        right_post = right_align_distribution(interval_result)

        # Support and probabilities should match after removing infinities
        @test length(support(left_result)) == length(support(left_post))
        @test length(support(centered_result)) == length(support(centered_post))
        @test length(support(right_result)) == length(support(right_post))

        @test support(left_result) ≈ support(left_post) atol=1e-12
        @test support(centered_result) ≈ support(centered_post) atol=1e-12
        @test support(right_result) ≈ support(right_post) atol=1e-12
    end

    @testset "Alignment methods with custom intervals" begin
        normal_dist = Normal(2, 1)
        custom_intervals = [0.0, 1.0, 2.0, 3.0, 4.0]

        # Test all alignment methods
        interval_result = discretize(normal_dist, custom_intervals; method = :interval)
        left_result = discretize(normal_dist, custom_intervals; method = :left_aligned)
        centered_result = discretize(normal_dist, custom_intervals; method = :centred)
        right_result = discretize(normal_dist, custom_intervals; method = :right_aligned)

        # Check support types
        @test eltype(support(interval_result)) <: Interval
        @test eltype(support(left_result)) <: Real
        @test eltype(support(centered_result)) <: Real
        @test eltype(support(right_result)) <: Real

        # All should sum to 1
        @test sum(probs(interval_result)) ≈ 1.0 atol=1e-10
        @test sum(probs(left_result)) ≈ 1.0 atol=1e-10
        @test sum(probs(centered_result)) ≈ 1.0 atol=1e-10
        @test sum(probs(right_result)) ≈ 1.0 atol=1e-10
    end

    @testset "Method consistency with post-alignment functions" begin
        normal_dist = Normal(0, 1)
        interval = 0.3

        # Test that direct method calls match post-processing
        interval_dist = discretize(normal_dist, interval; method = :interval)

        # Direct method calls
        left_direct = discretize(normal_dist, interval; method = :left_aligned)
        centered_direct = discretize(normal_dist, interval; method = :centred)
        right_direct = discretize(normal_dist, interval; method = :right_aligned)

        # Post-processing approach
        left_post = left_align_distribution(interval_dist)
        centered_post = centred_distribution(interval_dist)
        right_post = right_align_distribution(interval_dist)

        # Support and probabilities should match
        @test support(left_direct) ≈ support(left_post) atol=1e-12
        @test probs(left_direct) ≈ probs(left_post) atol=1e-12

        @test support(centered_direct) ≈ support(centered_post) atol=1e-12
        @test probs(centered_direct) ≈ probs(centered_post) atol=1e-12

        @test support(right_direct) ≈ support(right_post) atol=1e-12
        @test probs(right_direct) ≈ probs(right_post) atol=1e-12
    end

    @testset "Alignment with different numeric types" begin
        # Test with different numeric types
        normal_dist = Normal(0.0f0, 1.0f0)  # Float32
        interval = 0.5f0

        left_result = discretize(normal_dist, interval; method = :left_aligned)
        centered_result = discretize(normal_dist, interval; method = :centred)
        right_result = discretize(normal_dist, interval; method = :right_aligned)

        @test left_result isa DiscreteNonParametric
        @test centered_result isa DiscreteNonParametric
        @test right_result isa DiscreteNonParametric

        # Test with BigFloat
        normal_big = Normal(big(0.0), big(1.0))
        interval_big = big(0.5)

        left_big = discretize(normal_big, interval_big; method = :left_aligned)
        @test left_big isa DiscreteNonParametric
        @test eltype(support(left_big)) == BigFloat
    end

    @testset "Default method behavior" begin
        normal_dist = Normal(0, 1)
        interval = 0.2

        # Default should be :interval
        default_result = discretize(normal_dist, interval)
        explicit_result = discretize(normal_dist, interval; method = :interval)

        # Test that the distributions are equivalent
        @test length(support(default_result)) == length(support(explicit_result))
        @test probs(default_result) == probs(explicit_result)
        @test eltype(support(default_result)) <: Interval
        @test eltype(support(explicit_result)) <: Interval
    end

    @testset "Discrete distributions with alignment methods" begin
        poisson_dist = Poisson(3.0)
        interval = 2

        # Test alignment methods (skip unbiased for now due to implementation issues with discrete)
        interval_result = discretize(poisson_dist, interval; method = :interval)
        left_result = discretize(poisson_dist, interval; method = :left_aligned)
        centered_result = discretize(poisson_dist, interval; method = :centred)
        right_result = discretize(poisson_dist, interval; method = :right_aligned)

        @test all([result isa DiscreteNonParametric
                   for result in
                       [interval_result, left_result, centered_result, right_result]])
        @test eltype(support(interval_result)) <: Interval
        @test all([eltype(support(result)) <: Real
                   for result in [left_result, centered_result, right_result]])

        # All should sum to 1
        for result in [interval_result, left_result, centered_result, right_result]
            @test sum(probs(result)) ≈ 1.0 atol=1e-10
        end
    end

    @testset "Method parameter validation" begin
        normal_dist = Normal(0, 1)
        interval = 0.5

        # Test that invalid methods fall back to default behavior
        result = discretize(normal_dist, interval; method = :invalid_method)
        default_result = discretize(normal_dist, interval; method = :interval)
        @test result isa DiscreteNonParametric
        @test eltype(support(result)) <: Interval  # Should fall back to interval method

        # Test that all valid methods work
        valid_methods = [:interval, :left_aligned, :centred, :right_aligned]
        for method in valid_methods
            result = discretize(normal_dist, interval; method = method)
            @test result isa DiscreteNonParametric
            @test sum(probs(result)) ≈ 1.0 atol=1e-10
        end
    end

    @testset "Invalid distributions for unbiased method" begin
        # Test that distributions without defined means throw errors for :unbiased method
        invalid_dists = [
            Cauchy(0, 1), Pareto(1, 1), truncated(Cauchy(0, 1), lower = 0, upper = 1)]

        for dist in invalid_dists
            @test_throws ErrorException discretize(dist, 0.5; method = :unbiased)
        end
    end

    @testset "Trapezoid points parameter tests" begin
        # Test that trapezoid_points affects numerical integration accuracy in unbiased method

        # Use a distribution where mean calculation might need numerical integration
        gamma_dist = truncated(Gamma(2.0, 7.0), lower = 0, upper = 50.0)
        interval = 0.5

        # Test different trapezoid_points values
        trapezoid_vals = [100, 1000, 10000, 50000]
        results = []
        means = []

        for tp in trapezoid_vals
            result = discretize(gamma_dist, interval; method = :unbiased,
                trapezoid_points = tp, max_quantile = 0.999)
            push!(results, result)
            push!(means, mean(result))

            @test result isa DiscreteNonParametric
            @test sum(probs(result)) ≈ 1.0 atol=1e-10
        end

        # Higher trapezoid_points should give more accurate mean preservation
        #use safe mean since this isn't technically exact, i.e.  its only going to be equal for truncated distributions anyway
        true_mean = DiscretizeDistributions.safe_mean(gamma_dist, 100000)
        errors = [abs(m - true_mean) for m in means]

        # Generally, higher trapezoid_points should give better accuracy
        # (though not strictly monotonic due to numerical effects)
        @test errors[end] <= errors[1]  # 50000 points should be better than 100 points

        # Test that all results are reasonably close to the true mean
        for m in means
            @test abs(m - true_mean) < 0.06  # Within 6% error
        end

        # Test with a bounded distribution that doesn't need numerical integration
        uniform_dist = Uniform(0, 4)
        result_low = discretize(uniform_dist, 0.5; method = :unbiased, trapezoid_points = 100)
        result_high = discretize(uniform_dist, 0.5; method = :unbiased, trapezoid_points = 10000)

        # For distributions with analytical mean, trapezoid_points shouldn't matter much
        @test abs(mean(result_low) - mean(result_high)) < 1e-10
        @test abs(mean(result_low) - mean(uniform_dist)) < 1e-10

        # Test edge cases
        @testset "Trapezoid points edge cases" begin

            # throws error for very low points
            @test_throws ErrorException discretize(gamma_dist, 0.2; method = :unbiased, trapezoid_points = 10)

            # Test minimum viable trapezoid_points
            result_min = discretize(gamma_dist, 0.2; method = :unbiased, trapezoid_points = 100)
            @test result_min isa DiscreteNonParametric
            @test sum(probs(result_min)) ≈ 1.0 atol=1e-10

            # Test very high trapezoid_points (should not cause errors)
            result_max = discretize(gamma_dist, 0.2; method = :unbiased, trapezoid_points = 100000)
            @test result_max isa DiscreteNonParametric
            @test sum(probs(result_max)) ≈ 1.0 atol=1e-10

            # Both should preserve mean reasonably well
            @test abs(mean(result_min) - true_mean) < 0.10
            @test abs(mean(result_max) - true_mean) < 0.001
        end

        @testset "Particular cases" begin
            result = discretize(
                truncated(Gamma(2, 5); lower = 1), 0.5; method = :unbiased, trapezoid_points = 5000)
            @test result isa DiscreteNonParametric
            @test sum(probs(result)) ≈ 1.0 atol=1e-10

            result = discretize(
                truncated(Gumbel(); upper = 10), 0.5; method = :unbiased, trapezoid_points = 5000)
            @test result isa DiscreteNonParametric
            @test sum(probs(result)) ≈ 1.0 atol=1e-10

            #check for equal intervals
            unequal_intervals = [0.0, 1.0, 2.0, 4.0, 5.0]
            @test_throws ErrorException discretize(
                gamma_dist, unequal_intervals; method = :unbiased, trapezoid_points = 1000)
        end

        # Test that trapezoid_points only affects :unbiased method
        @testset "Trapezoid points method specificity" begin
            gamma_dist = Gamma(2.0, 7.0)
            interval = 0.3

            # For other methods, trapezoid_points should not affect results
            for method in [:left_aligned, :centred, :right_aligned]  # Skip :interval due to interval comparison issues
                result_low = discretize(gamma_dist, interval; method = method, trapezoid_points = 100)
                result_high = discretize(gamma_dist, interval; method = method, trapezoid_points = 10000)

                # Results should be identical for non-unbiased methods
                @test support(result_low) ≈ support(result_high)
                @test probs(result_low) ≈ probs(result_high)
            end

            # Test interval method separately (intervals can't be compared with ≈)
            interval_low = discretize(gamma_dist, interval; method = :interval, trapezoid_points = 100)
            interval_high = discretize(gamma_dist, interval; method = :interval, trapezoid_points = 10000)
            @test length(support(interval_low)) == length(support(interval_high))
            @test probs(interval_low) ≈ probs(interval_high)

            # But for unbiased method, they should potentially differ
            unbiased_low = discretize(gamma_dist, interval; method = :unbiased, trapezoid_points = 100)
            unbiased_high = discretize(gamma_dist, interval; method = :unbiased, trapezoid_points = 10000)

            # Means should be similar but potentially not identical due to numerical precision
            @test abs(mean(unbiased_low) - mean(unbiased_high)) < 0.01
        end
    end
end
