@doc """
    center_distribution(dist::Distributions.DiscreteNonParametric)

Center a discrete distribution by shifting support points to the midpoints between consecutive intervals.

This function takes the support points of a discrete distribution and creates a new distribution
where each support point is positioned at the center of the interval between consecutive original
support points. The last probability mass is dropped since there's no interval after the last point.

# Arguments
- `dist::Distributions.DiscreteNonParametric`: Input discrete distribution

# Returns
- `DiscreteNonParametric`: New distribution with centered support points

# Examples
```julia
using Distributions, DiscretiseDistributions

# Create a distribution with support at [1, 2, 3, 4]
dist = DiscreteNonParametric([1.0, 2.0, 3.0, 4.0], [0.25, 0.25, 0.25, 0.25])

# Center the distribution - support becomes [1.5, 2.5, 3.5]
centered = center_distribution(dist)
```
"""
function center_distribution(dist::Distributions.DiscreteNonParametric)
    if length(dist.support) <= 1
        error("Cannot center a distribution with one or no support points.")
    end
    xs = dist.support[1:(end - 1)] .+ (diff(dist.support) ./ 2)
    ps = dist.p[1:(end - 1)]
    return Distributions.DiscreteNonParametric(xs, ps; check_args=false)
end

@doc """
    center_distribution(dist::Distributions.DiscreteNonParametric, interval)

Center a discrete distribution by shifting all support points by half the given interval.

This function shifts all support points of the distribution by `interval / 2`, effectively
centering the distribution relative to the specified interval. All probability masses are preserved.

# Arguments
- `dist::Distributions.DiscreteNonParametric`: Input discrete distribution
- `interval`: The interval by which to center (support shifted by `interval / 2`)

# Returns
- `DiscreteNonParametric`: New distribution with shifted support points

# Examples
```julia
using Distributions, DiscretiseDistributions

# Create a distribution
dist = DiscreteNonParametric([1.0, 2.0, 3.0], [0.3, 0.4, 0.3])

# Center by interval of 1.0 - support becomes [1.5, 2.5, 3.5]
centered = center_distribution(dist, 1.0)
```
"""
function center_distribution(dist::Distributions.DiscreteNonParametric, interval)
    xs = dist.support .+ (interval / 2)
    return Distributions.DiscreteNonParametric(xs, dist.p)
end

@doc """
    right_align_distribution(dist::Distributions.DiscreteNonParametric)

Right-align a discrete distribution by shifting probabilities to the next support point.

This function creates a new distribution where each probability mass is associated with
the next support point (right-aligned). The last probability mass is dropped since
there's no support point to the right of the last one.

# Arguments
- `dist::Distributions.DiscreteNonParametric`: Input discrete distribution

# Returns
- `DiscreteNonParametric`: New distribution with right-aligned probabilities

# Examples
```julia
using Distributions, DiscretiseDistributions

# Create a distribution with support at [1, 2, 3, 4] and probs [0.2, 0.3, 0.3, 0.2]
dist = DiscreteNonParametric([1.0, 2.0, 3.0, 4.0], [0.2, 0.3, 0.3, 0.2])

# Right-align - support becomes [2, 3, 4] with probs [0.2, 0.3, 0.3]
right_aligned = right_align_distribution(dist)
```
"""
function right_align_distribution(dist::Distributions.DiscreteNonParametric)
    if length(dist.support) <= 1
        error("Cannot right-align a distribution with one or no support points.")
    end
    return Distributions.DiscreteNonParametric(dist.support[2:end], dist.p[1:(end - 1)]; check_args=false)
end

@doc """
    right_align_distribution(dist::Distributions.DiscreteNonParametric, interval)

Right-align a discrete distribution by shifting all support points by the given interval.

This function shifts all support points of the distribution by the specified interval,
effectively moving the entire distribution to the right. All probability masses are preserved.

# Arguments
- `dist::Distributions.DiscreteNonParametric`: Input discrete distribution
- `interval`: The interval by which to shift the support points

# Returns
- `DiscreteNonParametric`: New distribution with shifted support points

# Examples
```julia
using Distributions, DiscretiseDistributions

# Create a distribution
dist = DiscreteNonParametric([1.0, 2.0, 3.0], [0.3, 0.4, 0.3])

# Shift right by 0.5 - support becomes [1.5, 2.5, 3.5]
right_aligned = right_align_distribution(dist, 0.5)
```
"""
function right_align_distribution(dist::Distributions.DiscreteNonParametric, interval)
    xs = dist.support .+ (interval)
    return Distributions.DiscreteNonParametric(xs, dist.p)
end