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

"""
    discretise(dist::Distributions.DiscreteUnivariateDistribution, interval::Real; 
               min_quantile=0.001, max_quantile=0.999)

Discretise a discrete univariate distribution into intervals of fixed width.

This function groups the support of a discrete distribution into intervals of specified width,
aggregating probability masses within each interval. This is useful for reducing the granularity
of discrete distributions or for creating interval-based representations.

This is only sensible with discrete distributions that have finite support that makes sense as a continuous
distribution, such as Poisson or Binomial distributions. 

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

    return discretise(dist, xs)
end

"""
    discretise(dist::Distributions.DiscreteUnivariateDistribution, interval::AbstractVector)

Discretise a discrete univariate distribution using custom interval boundaries.

This function groups the support of a discrete distribution according to user-specified
interval boundaries, aggregating probability masses within each interval.

This is only sensible with discrete distributions that have finite support that makes sense as a continuous
distribution, such as Poisson or Binomial distributions. 


# Arguments
- `dist::Distributions.DiscreteUnivariateDistribution`: The discrete distribution to discretise
- `interval::AbstractVector`: Vector of interval boundaries (will be sorted automatically)

# Returns
- `DiscreteNonParametric`: Discrete distribution with aggregated probability masses

# Details
The input interval vector is automatically sorted. Probability masses are computed using the
probability density function (PDF) at floored interval boundaries. The distribution is assumed to be
uniform within the intervals of its original support.

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

    xs_floored = Int.(floor.(xs))

    #calculate full support and probability for the distribution over the range of xs
    full_support = collect(minimum(xs_floored):maximum(xs_floored))
    full_probability = map(Base.Fix1(Distributions.pdf, dist), full_support)
    full_probability /= sum(full_probability[1:(end-1)])

    n_intervals = length(xs) - 1

    probability = zeros(eltype(full_probability), n_intervals)

    for i in 1:n_intervals
        min_bound = xs[i]
        max_bound = xs[i + 1]
        #and in terms of the actual support of the distribution
        min_support = xs_floored[i]
        max_support = xs_floored[i + 1]

        #assume within that bound the distribution is uniform
        if min_support == max_support
            probability[i] = full_probability[min_support] * (max_bound - min_bound)
        else
            probability[i] =
                full_probability[min_support] * (min_support + 1 - min_bound) +
                full_probability[max_support] * (max_bound - max_support)
            if max_support > min_support + 1
                probability[i] += sum(full_probability[(min_support + 1):(max_support - 1)])
            end
        end
    end

    return Distributions.DiscreteNonParametric(xs[1:(end - 1)], probability)
end