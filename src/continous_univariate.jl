function define_bounds(dist::Distributions.ContinuousUnivariateDistribution, interval, min_quantile, max_quantile)
    max_dist_time = maximum(dist)
    addition = 0
    if isinf(max_dist_time)
        max_dist_time = Distributions.quantile(dist, max_quantile)
        #since its infinite it doesn't hurt to add an extra interval so we actually capture the quantile
        addition += 1
    end
    max_interval = Int((max_dist_time ÷ interval)) + addition

    min_dist_time = minimum(dist)
    subtraction = 0
    if isinf(min_dist_time)
        min_dist_time = Distributions.quantile(dist, min_quantile)
        subtraction += 1
    end
    min_interval = Int((min_dist_time ÷ interval)) - subtraction

    return min_interval:max_interval
end

function pseudo_cdf(dist::Distributions.ContinuousUnivariateDistribution, x::Real)
    return Distributions.cdf(dist, x - 1)
end

@doc """
    discretise(dist::Distributions.ContinuousUnivariateDistribution, interval::Real; 
               min_quantile=0.001, max_quantile=0.999)

Discretise a continuous univariate distribution into a discrete distribution using fixed intervals.

This function converts a continuous distribution into a discrete one by dividing the distribution's
support into intervals of fixed width and computing the probability mass in each interval using
the cumulative distribution function (CDF).

# Arguments
- `dist::Distributions.ContinuousUnivariateDistribution`: The continuous distribution to discretise
- `interval::Real`: The width of each discretisation interval
- `min_quantile=0.001`: Lower quantile bound for unbounded distributions
- `max_quantile=0.999`: Upper quantile bound for unbounded distributions

# Returns
- `DiscreteNonParametric`: Discrete distribution with probability masses at interval boundaries

# Details
For bounded distributions, the natural bounds are used. For unbounded distributions, the bounds
are determined using the specified quantiles. The probability mass in each interval is computed
as the difference in CDF values at the interval boundaries, ensuring the total probability sums to 1.

# Examples
```julia
using Distributions, DiscretiseDistributions

# Discretise a normal distribution with interval width 0.5
normal_dist = Normal(0, 1)
discrete_normal = discretise(normal_dist, 0.5)

# Discretise an exponential distribution with custom quantiles
exp_dist = Exponential(2.0)
discrete_exp = discretise(exp_dist, 0.1; min_quantile=0.01, max_quantile=0.95)
```
"""
function discretise(dist::Distributions.ContinuousUnivariateDistribution, interval::Real; min_quantile = 0.001, max_quantile=0.999)

    range = define_bounds(dist, interval, min_quantile, max_quantile)

    xs = collect(range * interval)

    return discretise(dist, xs)
end

@doc """
    discretise(dist::Distributions.UnivariateDistribution, interval::AbstractVector)

Discretise a univariate distribution using custom interval boundaries.

This function converts a continuous distribution into a discrete one using user-specified
interval boundaries. The probability mass in each interval is computed using the cumulative
distribution function (CDF) or pseudo CDF (which assume a discrete distribution only supports integers
and is uniform between support points).

# Arguments
- `dist::Distributions.UnivariateDistribution`: The distribution to discretise
- `interval::AbstractVector`: Vector of interval boundaries (will be sorted automatically)

# Returns
- `DiscreteNonParametric`: Discrete distribution with probability masses at specified boundaries

# Details
The input interval vector is automatically sorted. Probability masses are computed as differences
in CDF values between consecutive boundaries. The resulting discrete distribution has support
points at the interval boundaries (excluding the last boundary) with corresponding probabilities.

# Examples
```julia
using Distributions, DiscretiseDistributions

# Discretise using custom intervals
normal_dist = Normal(5, 2)
custom_intervals = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
discrete_normal = discretise(normal_dist, custom_intervals)

# The intervals will be sorted automatically if needed
unsorted_intervals = [8.0, 0.0, 4.0, 2.0, 10.0]
discrete_normal2 = discretise(normal_dist, unsorted_intervals)

# Discrete distribution
poisson_dist = Poisson(3.0)
discrete_poisson = discretise(poisson_dist, [0.5, 2, 4, 6, 8, 10])
```
"""
function discretise(dist::Distributions.UnivariateDistribution, interval::AbstractVector)

    xs = sort(interval)

    probability = diff(pseudo_cdf.(dist, xs))

    xs = xs[1:(end-1)]

    #set to sum to one
    probability /= sum(probability)

    return Distributions.DiscreteNonParametric(xs, probability)
end