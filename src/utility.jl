function remove_infinities(
        dist::Distributions.DiscreteNonParametric{
            IntervalArithmetic.Interval{T}, P, S, V
        }) where {T <: Real, P <: Real, S <: AbstractVector{IntervalArithmetic.Interval{T}}, V <: AbstractVector{P}}

    support = dist.support
    p = dist.p
    if isinf(IntervalArithmetic.inf(support[1]))
        @warn "Support contains interval with negative infinity, removing."

        support = support[2:end]
        p = p[2:end]
    end
    if isinf(IntervalArithmetic.sup(support[end]))
        @warn "Support contains interval with positive infinity, removing."

        support = support[1:(end - 1)]
        p = p[1:(end - 1)]
    end
    p = p ./ sum(p)  # Normalize probabilities
    return Distributions.DiscreteNonParametric(support, p; check_args=false)
end

@doc """
    centred_distribution(dist::Distributions.DiscreteNonParametric{IntervalArithmetic.Interval{T}, ...})

Convert an interval-based discrete distribution to a centered point-based distribution.

This function takes a discrete distribution with interval support and creates a new distribution
where each support point is positioned at the center (midpoint) of the corresponding interval.
Infinite intervals are automatically removed before conversion.

# Arguments
- `dist::DiscreteNonParametric{Interval{T}, ...}`: Input discrete distribution with interval support

# Returns
- `DiscreteNonParametric{T, ...}`: New distribution with centered point support

# Examples
```julia
using Distributions, DiscretizeDistributions, IntervalArithmetic

# Create an interval-based distribution
intervals = [interval(0.0, 1.0), interval(1.0, 2.0), interval(2.0, 3.0)]
probs = [0.3, 0.4, 0.3]
interval_dist = DiscreteNonParametric(intervals, probs, check_args=false)

# Convert to centered points
centered = centred_distribution(interval_dist)
# Support becomes [0.5, 1.5, 2.5] (midpoints of intervals)
```
"""
function centred_distribution(
    dist::Distributions.DiscreteNonParametric{
        IntervalArithmetic.Interval{T}, P, S, V
    }) where {T <: Real, P <: Real, S <: AbstractVector{IntervalArithmetic.Interval{T}}, V <: AbstractVector{P}}
    
    dist = remove_infinities(dist)

    xs = (IntervalArithmetic.sup.(dist.support) .+ IntervalArithmetic.inf.(dist.support)) ./ 2

    return Distributions.DiscreteNonParametric(xs, dist.p)
end

@doc """
    left_align_distribution(dist::Distributions.DiscreteNonParametric{IntervalArithmetic.Interval{T}, ...})

Convert an interval-based discrete distribution to a left-aligned point-based distribution.

This function takes a discrete distribution with interval support and creates a new distribution
where each support point is positioned at the left endpoint (infimum) of the corresponding interval.
Infinite intervals are automatically removed before conversion.

# Arguments
- `dist::DiscreteNonParametric{Interval{T}, ...}`: Input discrete distribution with interval support

# Returns
- `DiscreteNonParametric{T, ...}`: New distribution with left-aligned point support

# Examples
```julia
using Distributions, DiscretizeDistributions, IntervalArithmetic

# Create an interval-based distribution
intervals = [interval(0.0, 1.0), interval(1.0, 2.0), interval(2.0, 3.0)]
probs = [0.3, 0.4, 0.3]
interval_dist = DiscreteNonParametric(intervals, probs, check_args=false)

# Convert to left-aligned points
left_aligned = left_align_distribution(interval_dist)
# Support becomes [0.0, 1.0, 2.0] (left endpoints of intervals)
```
"""
function left_align_distribution(
    dist::Distributions.DiscreteNonParametric{
        IntervalArithmetic.Interval{T}, P, S, V
    }) where {T <: Real, P <: Real, S <: AbstractVector{IntervalArithmetic.Interval{T}}, V <: AbstractVector{P}}
    
    dist = remove_infinities(dist)

    xs = IntervalArithmetic.inf.(dist.support)

    return Distributions.DiscreteNonParametric(xs, dist.p)
end

@doc """
    right_align_distribution(dist::Distributions.DiscreteNonParametric{IntervalArithmetic.Interval{T}, ...})

Convert an interval-based discrete distribution to a right-aligned point-based distribution.

This function takes a discrete distribution with interval support and creates a new distribution
where each support point is positioned at the right endpoint (supremum) of the corresponding interval.
Infinite intervals are automatically removed before conversion.

# Arguments
- `dist::DiscreteNonParametric{Interval{T}, ...}`: Input discrete distribution with interval support

# Returns
- `DiscreteNonParametric{T, ...}`: New distribution with right-aligned point support

# Examples
```julia
using Distributions, DiscretizeDistributions, IntervalArithmetic

# Create an interval-based distribution
intervals = [interval(0.0, 1.0), interval(1.0, 2.0), interval(2.0, 3.0)]
probs = [0.3, 0.4, 0.3]
interval_dist = DiscreteNonParametric(intervals, probs, check_args=false)

# Convert to right-aligned points
right_aligned = right_align_distribution(interval_dist)
# Support becomes [1.0, 2.0, 3.0] (right endpoints of intervals)
```
"""
function right_align_distribution(dist::Distributions.DiscreteNonParametric{
        IntervalArithmetic.Interval{T}, P, S, V
    }) where {T <: Real, P <: Real, S <: AbstractVector{IntervalArithmetic.Interval{T}}, V <: AbstractVector{P}}
    
    dist = remove_infinities(dist)

    xs = IntervalArithmetic.sup.(dist.support)

    return Distributions.DiscreteNonParametric(xs, dist.p)
end