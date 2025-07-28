function define_bounds(dist::Distributions.DiscreteUnivariateDistribution, interval, min_quantile, max_quantile)
    max_dist_time = maximum(dist)
    if isinf(max_dist_time)
        max_dist_time = Distributions.quantile(dist, max_quantile) + 1
    end
    max_interval = Int((max_dist_time ÷ interval))

    min_dist_time = minimum(dist)
    if isinf(min_dist_time)
        min_dist_time = Distributions.quantile(dist, min_quantile) + 1
    end
    min_interval = Int((min_dist_time ÷ interval))

    return min_interval:max_interval
end

function discretise_univariate_discrete(dist, xs)
    xs_daily = floor.(xs)
    probability = map(Base.Fix1(Distributions.pdf, dist), xs_daily)

    #set to sum to one
    probability /= sum(probability)

    return Distributions.DiscreteNonParametric(xs[probability .> 0], probability[probability .> 0])
end

@doc """
    discretise(dist::Distributions.DiscreteUnivariateDistribution, interval::Real; 
               min_quantile=0.001, max_quantile=0.999)

Discretise a discrete univariate distribution into intervals of fixed width.

This function groups the support of a discrete distribution into intervals of specified width,
aggregating probability masses within each interval. This is useful for reducing the granularity
of discrete distributions or for creating interval-based representations.

# Arguments
- `dist::Distributions.DiscreteUnivariateDistribution`: The discrete distribution to discretise
- `interval::Real`: The width of each discretisation interval
- `min_quantile=0.001`: Lower quantile bound for unbounded distributions
- `max_quantile=0.999`: Upper quantile bound for unbounded distributions

# Returns
- `DiscreteNonParametric`: Discrete distribution with aggregated probability masses

# Details
For bounded distributions, the natural bounds are used. For unbounded distributions, bounds
are determined using the specified quantiles. Probability masses are computed using the
probability density function (PDF) at floored interval boundaries and normalized to sum to 1.

# Examples
```julia
using Distributions, DiscretiseDistributions

# Discretise a Poisson distribution with interval width 2
poisson_dist = Poisson(5.0)
discrete_poisson = discretise(poisson_dist, 2)

# Discretise a Binomial distribution
binomial_dist = Binomial(20, 0.3)
discrete_binomial = discretise(binomial_dist, 3; max_quantile=0.95)
```
"""
function discretise(dist::Distributions.DiscreteUnivariateDistribution, interval::Real; min_quantile = 0.001, max_quantile=0.999)

    range = define_bounds(dist, interval, min_quantile, max_quantile)
    
    xs = collect(range * interval)

    return discretise_univariate_discrete(dist, xs)
end

@doc """
    discretise(dist::Distributions.DiscreteUnivariateDistribution, interval::AbstractVector)

Discretise a discrete univariate distribution using custom interval boundaries.

This function groups the support of a discrete distribution according to user-specified
interval boundaries, aggregating probability masses within each interval.

# Arguments
- `dist::Distributions.DiscreteUnivariateDistribution`: The discrete distribution to discretise
- `interval::AbstractVector`: Vector of interval boundaries (will be sorted automatically)

# Returns
- `DiscreteNonParametric`: Discrete distribution with aggregated probability masses

# Details
The input interval vector is automatically sorted. Probability masses are computed using the
probability density function (PDF) at floored interval boundaries. The resulting distribution
has support points at the interval boundaries with aggregated probabilities normalized to sum to 1.

# Examples
```julia
using Distributions, DiscretiseDistributions

# Discretise using custom intervals
poisson_dist = Poisson(3.0)
custom_intervals = [0, 2, 4, 6, 10, 15]
discrete_poisson = discretise(poisson_dist, custom_intervals)

# Intervals are sorted automatically
unsorted_intervals = [10, 0, 4, 2, 15]
discrete_poisson2 = discretise(poisson_dist, unsorted_intervals)
```
"""
function discretise(dist::Distributions.DiscreteUnivariateDistribution, interval::AbstractVector)

    xs = sort(interval)

    return discretise_univariate_discrete(dist, xs)
end