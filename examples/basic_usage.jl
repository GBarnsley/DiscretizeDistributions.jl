#!/usr/bin/env julia

"""
Example script demonstrating DiscretiseDistributions.jl functionality

This script shows how to use the package to discretise various types of distributions
and apply different alignment transformations.
"""

using Distributions
using DiscretiseDistributions

println("=== DiscretiseDistributions.jl Examples ===\n")

# Example 1: Discretising a Normal distribution
println("1. Discretising a Normal distribution")
println("   Original: Normal(μ=0, σ=1)")

normal_dist = Normal(0, 1)
discrete_normal = discretise(normal_dist, 0.5)

println("   Discretised with interval 0.5:")
println("   Support: ", support(discrete_normal))
println("   Probabilities: ", round.(probs(discrete_normal), digits=4))
println("   Total probability: ", sum(probs(discrete_normal)))
println()

# Example 2: Using custom intervals
println("2. Using custom intervals")
custom_intervals = [-3.0, -1.0, -0.5, 0.0, 0.5, 1.0, 3.0]
discrete_custom = discretise(normal_dist, custom_intervals)

println("   Custom intervals: ", custom_intervals)
println("   Support: ", support(discrete_custom))
println("   Probabilities: ", round.(probs(discrete_custom), digits=4))
println()

# Example 3: Discretising an Exponential distribution
println("3. Discretising an Exponential distribution")
println("   Original: Exponential(λ=2.0)")

exp_dist = Exponential(2.0)
discrete_exp = discretise(exp_dist, 0.2; max_quantile=0.95)

println("   Discretised with interval 0.2 (95% quantile bound):")
println("   Support: ", support(discrete_exp)[1:min(10, end)], length(support(discrete_exp)) > 10 ? "..." : "")
println("   First 5 probabilities: ", round.(probs(discrete_exp)[1:min(5, end)], digits=4))
println()

# Example 4: Discretising a discrete distribution (Poisson)
println("4. Discretising a Poisson distribution")
println("   Original: Poisson(λ=3.0)")

poisson_dist = Poisson(3.0)
discrete_poisson = discretise(poisson_dist, 2)  # Group into intervals of width 2

println("   Grouped into intervals of width 2:")
println("   Support: ", support(discrete_poisson))
println("   Probabilities: ", round.(probs(discrete_poisson), digits=4))
println()

# Example 5: Distribution alignment transformations
println("5. Distribution alignment transformations")

# Create a simple discrete distribution
test_dist = DiscreteNonParametric([1.0, 2.0, 3.0, 4.0], [0.2, 0.3, 0.3, 0.2])
println("   Original distribution:")
println("   Support: ", support(test_dist))
println("   Probabilities: ", probs(test_dist))

# Center the distribution
centered = center_distribution(test_dist)
println("   Centered distribution:")
println("   Support: ", support(centered))
println("   Probabilities: ", probs(centered))

# Right-align the distribution
right_aligned = right_align_distribution(test_dist)
println("   Right-aligned distribution:")
println("   Support: ", support(right_aligned))
println("   Probabilities: ", probs(right_aligned))

# Shift by custom amount
shifted = right_align_distribution(test_dist, 0.5)
println("   Shifted by 0.5:")
println("   Support: ", support(shifted))
println("   Probabilities: ", probs(shifted))
println()

# Example 6: Working with bounded distributions
println("6. Working with bounded distributions")
println("   Original: Uniform(0, 10)")

uniform_dist = Uniform(0, 10)
discrete_uniform = discretise(uniform_dist, 1.0)

println("   Discretised with interval 1.0:")
println("   Support: ", support(discrete_uniform))
println("   Probabilities: ", round.(probs(discrete_uniform), digits=4))
println("   Note: All probabilities are equal for uniform distribution")
println()

println("=== Examples completed ===")
