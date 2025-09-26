function define_intervals(
        dist::Distributions.UnivariateDistribution, interval, min_quantile, max_quantile)
    #to avoid issues with precision we'll translate everything into multiples of the interval

    min_value = minimum(dist)
    if isinf(min_value)
        #set next value to the interval of given size that contains the minimum quantile
        finite_minimum = Int((Distributions.quantile(dist, min_quantile) ÷ interval)) - 1
    else
        finite_minimum = Int(min_value ÷ interval)
    end

    max_value = maximum(dist)
    if isinf(max_value)
        #set finite max to the interval of given size that contains the maximum quantile
        finite_maximum = Int((Distributions.quantile(dist, max_quantile) ÷ interval)) + 1
    else
        finite_maximum = Int(max_value ÷ interval)
    end

    lower_bound = interval .* collect(finite_minimum:finite_maximum) #ensure we have a zero in the lower bound

    #final additions to capture the tails
    if min_value != lower_bound[1]
        lower_bound = vcat([min_value], lower_bound)
    end

    if max_value != lower_bound[end]
        upper_bound = vcat(lower_bound[2:end], [max_value])
    else
        upper_bound = lower_bound[2:end]
        lower_bound = lower_bound[1:(end - 1)]
    end

    return IntervalArithmetic.interval.(lower_bound, upper_bound)
end

function safe_mean(dist::Distributions.UnivariateDistribution, trapezoid_points)
    try
        return Distributions.mean(dist)
    catch
        # Fallback to numerical integration: E[X] = ∫ x * f(x) dx
        min_val = minimum(dist)
        max_val = maximum(dist)

        # Handle infinite bounds
        if isinf(min_val)
            min_val = Distributions.quantile(dist, 1e-12)
        end
        if isinf(max_val)
            max_val = Distributions.quantile(dist, 1 - 1e-12)
        end

        # Simple numerical integration using trapezoidal rule
        x_vals = range(min_val, max_val; length = trapezoid_points)

        integrand_vals = x_vals .* Distributions.pdf.(Ref(dist), x_vals)
        dx = (max_val - min_val) / (trapezoid_points - 1)

        # Trapezoidal rule
        return dx * (0.5 * (integrand_vals[1] + integrand_vals[end]) + sum(integrand_vals[2:end-1]))
    end
end

function limited_expectation(dist::Distributions.UnivariateDistribution, limit::Real, trapezoid_points)
    # Handle edge cases
    if limit >= maximum(dist)
        # If limit is at or above the maximum possible value, return the mean
        return safe_mean(dist, trapezoid_points)
    elseif limit <= minimum(dist)
        # If limit is at or below minimum, return the limit
        return limit
    else
        # For the general case: E[min(X, u)] = E[X·1_{X≤u}] + u·P(X > u)
        # This equals: ∫₋∞ᵘ x·f(x)dx + u·P(X > u)

        # Calculate P(X > u)
        survival_prob = Distributions.ccdf(dist, limit)

        # Calculate E[X·1_{X≤u}] using truncated distribution
        if survival_prob ≈ 1.0
            # Essentially all mass is above the limit
            return limit
        elseif survival_prob ≈ 0.0
            # Essentially all mass is below the limit
            return safe_mean(dist, trapezoid_points)
        else
            truncated_dist = Distributions.truncated(dist; lower=minimum(dist), upper=limit)
            truncated_expectation = safe_mean(truncated_dist, trapezoid_points)
            truncated_prob = 1.0 - survival_prob
            return truncated_expectation * truncated_prob + limit * survival_prob
        end
    end
end

function pseudo_cdf(dist::Distributions.ContinuousUnivariateDistribution, x::Real)
    return Distributions.cdf(dist, x)
end

function pseudo_cdf(dist::Distributions.DiscreteUnivariateDistribution, x::Real)
    # Handle infinite values
    if isinf(x) && x > 0
        return 1.0
    elseif isinf(x) && x < 0
        return 0.0
    end

    floor_x = floor(x)
    return Distributions.cdf(dist, floor_x - 1) +
           Distributions.pdf(dist, floor_x) * (x - floor_x)
end

@doc """
    discretize(dist::Distributions.UnivariateDistribution, interval::Real;
               method=:interval, min_quantile=0.001, max_quantile=0.999)

Discretize a univariate distribution into a discrete distribution using fixed intervals.

This function converts a univariate distribution into a discrete one by dividing the distribution's
support into intervals of fixed width and computing the probability mass in each interval.

# Arguments
- `dist::Distributions.UnivariateDistribution`: The distribution to discretize (continuous or discrete)
- `interval::Real`: The width of each discretisation interval
- `method::Symbol=:interval`: Method for representing the output distribution
    - `:interval` (default): Return `IntervalArithmetic.Interval` objects as support
    - `:left_aligned`: Convert intervals to left-aligned point values
    - `:centred`: Convert intervals to centered point values
    - `:right_aligned`: Convert intervals to right-aligned point values
    - `:unbiased`: Return unbiased point estimates (requires equal interval widths), designed such that the means match, see `discretize` from the R package `actuar`
- `min_quantile=0.001`: Lower quantile bound for unbounded distributions
- `max_quantile=0.999`: Upper quantile bound for unbounded distributions
- `trapezoid_points::Int=10000`: Number of points for numerical integration of the mean (when needed)

# Returns
- `DiscreteNonParametric`: Discrete distribution with support determined by the method parameter

# Details
For bounded distributions, the natural bounds are used. For unbounded distributions, the bounds
are determined using the specified quantiles. The probability mass in each interval is computed
using the CDF for continuous distributions or a pseudo-CDF for discrete distributions.

# Examples
```julia
using Distributions, DiscretizeDistributions, IntervalArithmetic

# Discretize a normal distribution with interval width 0.5
normal_dist = Normal(0, 1)

# Different output methods
discrete_intervals = discretize(normal_dist, 0.5)                        # Intervals (default)
discrete_left = discretize(normal_dist, 0.5; method=:left_aligned)      # Left endpoints
discrete_center = discretize(normal_dist, 0.5; method=:centred)         # Midpoints
discrete_right = discretize(normal_dist, 0.5; method=:right_aligned)    # Right endpoints

# Compare means (centered method typically closest to original)
println("Original mean: ", mean(normal_dist))
println("Centered discretization mean: ", mean(discrete_center))

# Discretize a discrete distribution
poisson_dist = Poisson(3.0)
discrete_poisson = discretize(poisson_dist, 2; method=:centred)
```
"""
dist = Distributions.Normal(0.0, 1.0); interval = 1//10
function discretize(dist::Distributions.UnivariateDistribution, interval::Real; method::Symbol = :interval,
        min_quantile = 0.001, max_quantile = 0.999, trapezoid_points::Int = 10000)
    xs = define_intervals(dist, interval, min_quantile, max_quantile)

    return discretize(dist, xs; method = method, trapezoid_points = trapezoid_points)
end

@doc """
    discretize(dist::Distributions.UnivariateDistribution, interval::AbstractVector; method=:interval)

Discretize a univariate distribution using custom interval boundaries.

This function converts a univariate distribution into a discrete one using user-specified
interval boundaries. The resulting distribution represents the probability mass in each interval.

# Arguments
- `dist::Distributions.UnivariateDistribution`: The distribution to discretize
- `interval::AbstractVector`: Vector of interval boundaries (will be sorted automatically)
- `method::Symbol=:interval`: Method for representing the output distribution
    - `:interval` (default): Return `IntervalArithmetic.Interval` objects as support
    - `:left_aligned`: Convert intervals to left-aligned point values
    - `:centred`: Convert intervals to centered point values
    - `:right_aligned`: Convert intervals to right-aligned point values
    - `:unbiased`: Return unbiased point estimates (requires equal interval widths), designed such that the means match, see `discretize` from the R package `actuar`
- `trapezoid_points::Int=10000`: Number of points for numerical integration of the mean (when needed)

# Returns
- `DiscreteNonParametric`: Discrete distribution with support determined by the method parameter

# Details
The input interval vector is automatically sorted and combined with distribution bounds.
Probability masses are computed using the CDF for continuous distributions or pseudo-CDF
for discrete distributions. The resulting distribution represents probability masses over
intervals `[a_i, a_{i+1})`.

For the `:unbiased` method with unequal intervals, the function will warn and fall back to `:centred`.

# Examples
```julia
using Distributions, DiscretizeDistributions, IntervalArithmetic

# Discretize using custom intervals
normal_dist = Normal(5, 2)
custom_intervals = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]

# Different output methods
discrete_intervals = discretize(normal_dist, custom_intervals)                      # Intervals
discrete_left = discretize(normal_dist, custom_intervals; method=:left_aligned)    # Left points
discrete_center = discretize(normal_dist, custom_intervals; method=:centred)       # Midpoints
discrete_right = discretize(normal_dist, custom_intervals; method=:right_aligned)  # Right points

# Support: intervals like [interval(-∞, 0.0), interval(0.0, 2.0), ..., interval(10.0, ∞)]

# Discrete distribution with custom intervals
poisson_dist = Poisson(3.0)
discrete_poisson = discretize(poisson_dist, [0.5, 2, 4, 6, 8, 10]; method=:centred)
```
"""
function discretize(dist::Distributions.UnivariateDistribution, interval::AbstractVector; method::Symbol = :interval, trapezoid_points::Int = 10000)
    xs = sort(interval)
    if method == :unbiased
        # Check if intervals are equally spaced
        if length(xs) >= 2
            interval_widths = diff(xs)
            if length(interval_widths) >= 2 && any([!(i ≈ interval_widths[1]) for i in interval_widths[2:end]])
                @error "The :unbiased method is not compatible with unequal interval sizes."
            end
        end
    end

    min_value = minimum(dist)

    max_value = maximum(dist)

    xs = xs[xs .> min_value]

    xs = xs[xs .< max_value]

    xs = vcat([min_value], xs, [max_value])

    xs = IntervalArithmetic.interval.(xs[1:(end - 1)], xs[2:end])

    return discretize(dist, xs; method = method, trapezoid_points = trapezoid_points)
end

@doc """
    discretize(dist::Distributions.UnivariateDistribution,
               interval::AbstractVector{IntervalArithmetic.Interval{X}}; method=:interval) where X <: Real

Discretize a univariate distribution using pre-constructed interval objects.

This function converts a univariate distribution into a discrete one using user-specified
`IntervalArithmetic.Interval` objects. This is the core discretization method that all other
`discretize` methods ultimately call.

# Arguments
- `dist::Distributions.UnivariateDistribution`: The distribution to discretize
- `interval::AbstractVector{IntervalArithmetic.Interval{X}}`: Vector of pre-constructed intervals
- `method::Symbol=:interval`: Method for representing the output distribution
    - `:interval` (default): Return intervals as support points
    - `:left_aligned`: Convert intervals to left-aligned point values
    - `:centred`: Convert intervals to centered point values
    - `:right_aligned`: Convert intervals to right-aligned point values
    - `:unbiased`: Return unbiased point estimates (requires equal interval widths), designed such that the means match, see `discretize` from the R package `actuar`
- `trapezoid_points::Int=10000`: Number of points for numerical integration of the mean (when needed)

# Returns
- `DiscreteNonParametric`: Discrete distribution with support determined by the method parameter

# Details
This method computes probability masses directly using the interval boundaries. For each interval
`[a, b]`, the probability is computed as `cdf(dist, b) - cdf(dist, a)`. The
resulting probabilities are normalized to sum to 1.

# Examples
```julia
using Distributions, DiscretizeDistributions, IntervalArithmetic

# Create intervals manually
intervals = [interval(-1.0, 0.0), interval(0.0, 1.0), interval(1.0, 2.0)]

# Discretize using these intervals with different methods
normal_dist = Normal(0, 1)
discrete_intervals = discretize(normal_dist, intervals)                      # Intervals
discrete_centered = discretize(normal_dist, intervals; method=:centred)      # Midpoints
discrete_left = discretize(normal_dist, intervals; method=:left_aligned)     # Left endpoints

# Each method gives the same probabilities but different support representations
```
"""
function discretize(dist::Distributions.UnivariateDistribution,
        interval::AbstractVector{IntervalArithmetic.Interval{X}}; method::Symbol = :interval, trapezoid_points::Int = 10000) where {X <: Real}
    xs_values = vcat(
        IntervalArithmetic.inf.(interval), [IntervalArithmetic.sup(interval[end])])

    if method == :unbiased

        unbiased_values = xs_values[.!isinf.(xs_values)]

        interval_widths = diff(unbiased_values)
        if length(interval_widths) >= 2 && any([!(i ≈ interval_widths[1]) for i in interval_widths[2:end]])
            error("The :unbiased method requires equal interval widths.")
        end

        probability = Vector{eltype(dist)}(undef, length(unbiased_values))

        probability[1] = (limited_expectation(dist, unbiased_values[1], trapezoid_points) - limited_expectation(dist, unbiased_values[2], trapezoid_points))/(unbiased_values[2] - unbiased_values[1]) + Distributions.ccdf(dist, unbiased_values[1])
        probability[end] = (limited_expectation(dist, unbiased_values[end], trapezoid_points) - limited_expectation(dist, unbiased_values[end - 1], trapezoid_points))/(unbiased_values[end] - unbiased_values[end - 1]) - Distributions.ccdf(dist, unbiased_values[end])

        probability[2:(end - 1)] .= ((2 .* limited_expectation.(dist, unbiased_values[2:(end - 1)], trapezoid_points)) .- limited_expectation.(dist, unbiased_values[1:(end - 2)], trapezoid_points) .- limited_expectation.(dist, unbiased_values[3:end], trapezoid_points))./
            (unbiased_values[3:end] .- unbiased_values[2:(end - 1)])
    else
        probability = diff(pseudo_cdf.(dist, xs_values))
    end

    #set to sum to one
    probability /= sum(probability)

    if method == :right_aligned
        return right_align_distribution(Distributions.DiscreteNonParametric(interval, probability; check_args = false))
    elseif method == :left_aligned
        return left_align_distribution(Distributions.DiscreteNonParametric(interval, probability; check_args = false))
    elseif method == :centred
        return centred_distribution(Distributions.DiscreteNonParametric(interval, probability; check_args = false))
    elseif method == :unbiased
        return Distributions.DiscreteNonParametric(unbiased_values, probability)
    else
        return Distributions.DiscreteNonParametric(interval, probability; check_args = false)
    end
end
