function define_intervals(dist::Distributions.UnivariateDistribution, interval, min_quantile, max_quantile)
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
        lower_bound = lower_bound[1:(end-1)]
    end

    return IntervalArithmetic.interval.(lower_bound, upper_bound)
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
    return Distributions.cdf(dist, floor_x - 1) + Distributions.pdf(dist, floor_x) * (x - floor_x)
end

@doc """
    discretise(dist::Distributions.UnivariateDistribution, interval::Real; 
               min_quantile=0.001, max_quantile=0.999)

Discretise a univariate distribution into an interval-based discrete distribution using fixed intervals.

This function converts a univariate distribution into a discrete one by dividing the distribution's
support into intervals of fixed width and computing the probability mass in each interval. The 
resulting distribution has `IntervalArithmetic.Interval` objects as support points.

# Arguments
- `dist::Distributions.UnivariateDistribution`: The distribution to discretise (continuous or discrete)
- `interval::Real`: The width of each discretisation interval
- `min_quantile=0.001`: Lower quantile bound for unbounded distributions
- `max_quantile=0.999`: Upper quantile bound for unbounded distributions

# Returns
- `DiscreteNonParametric{Interval{T}, ...}`: Discrete distribution with interval support

# Details
For bounded distributions, the natural bounds are used. For unbounded distributions, the bounds
are determined using the specified quantiles. The probability mass in each interval is computed
using the CDF for continuous distributions or a pseudo-CDF for discrete distributions.

# Examples
```julia
using Distributions, DiscretiseDistributions, IntervalArithmetic

# Discretise a normal distribution with interval width 0.5
normal_dist = Normal(0, 1)
discrete_intervals = discretise(normal_dist, 0.5)
# Returns distribution with support like [interval(-∞, -3.5), interval(-3.5, -3.0), ...]

# Convert to point-based distributions if needed
left_aligned = left_align_distribution(discrete_intervals)    # [-3.5, -3.0, -2.5, ...]
centered = centred_distribution(discrete_intervals)          # [-3.25, -2.75, -2.25, ...]

# Discretise a discrete distribution
poisson_dist = Poisson(3.0)
discrete_poisson = discretise(poisson_dist, 2)
```
"""
function discretise(dist::Distributions.UnivariateDistribution, interval::Real; min_quantile = 0.001, max_quantile=0.999)

    xs = define_intervals(dist, interval, min_quantile, max_quantile)

    return discretise(dist, xs)
end

@doc """
    discretise(dist::Distributions.UnivariateDistribution, interval::AbstractVector)

Discretise a univariate distribution using custom interval boundaries.

This function converts a univariate distribution into a discrete one using user-specified
interval boundaries. The resulting distribution has `IntervalArithmetic.Interval` objects 
as support points representing the probability mass in each interval.

# Arguments
- `dist::Distributions.UnivariateDistribution`: The distribution to discretise
- `interval::AbstractVector`: Vector of interval boundaries (will be sorted automatically)

# Returns
- `DiscreteNonParametric{Interval{T}, ...}`: Discrete distribution with interval support

# Details
The input interval vector is automatically sorted and combined with distribution bounds.
Probability masses are computed using the CDF for continuous distributions or pseudo-CDF
for discrete distributions. The resulting distribution represents probability masses over
intervals `[a_i, a_{i+1})`.

# Examples
```julia
using Distributions, DiscretiseDistributions, IntervalArithmetic

# Discretise using custom intervals
normal_dist = Normal(5, 2)
custom_intervals = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
discrete_intervals = discretise(normal_dist, custom_intervals)
# Support: [interval(-∞, 0.0), interval(0.0, 2.0), ..., interval(10.0, ∞)]

# Convert to different alignments
left_points = left_align_distribution(discrete_intervals)
centered_points = centred_distribution(discrete_intervals)
right_points = right_align_distribution(discrete_intervals)

# Discrete distribution with custom intervals
poisson_dist = Poisson(3.0)
discrete_poisson = discretise(poisson_dist, [0.5, 2, 4, 6, 8, 10])
```
"""
function discretise(dist::Distributions.UnivariateDistribution, interval::AbstractVector)
    
    xs = sort(interval)

    min_value = minimum(dist)

    max_value = maximum(dist)

    xs = xs[xs .> min_value]

    xs = xs[xs .< max_value]

    xs = vcat([min_value], xs, [max_value])

    xs = IntervalArithmetic.interval.(xs[1:(end-1)], xs[2:end])

    return discretise(dist, xs)
end

@doc """
    discretise(dist::Distributions.UnivariateDistribution, interval::AbstractVector{IntervalArithmetic.Interval{X}}) where X <: Real

Discretise a univariate distribution using pre-constructed interval objects.

This function converts a univariate distribution into a discrete one using user-specified
`IntervalArithmetic.Interval` objects. This is the core discretization method that all other
`discretise` methods ultimately call. The resulting distribution has the same interval objects
as support points with computed probability masses.

# Arguments
- `dist::Distributions.UnivariateDistribution`: The distribution to discretise
- `interval::AbstractVector{IntervalArithmetic.Interval{X}}`: Vector of pre-constructed intervals

# Returns
- `DiscreteNonParametric{Interval{X}, ...}`: Discrete distribution with interval support

# Details
This method computes probability masses directly using the interval boundaries. For each interval
`[a, b]`, the probability is computed as `cdf(dist, b) - cdf(dist, a)`. The 
resulting probabilities are normalized to sum to 1.

# Examples
```julia
using Distributions, DiscretiseDistributions, IntervalArithmetic

# Create intervals manually
intervals = [interval(-1.0, 0.0), interval(0.0, 1.0), interval(1.0, 2.0)]

# Discretise using these intervals
normal_dist = Normal(0, 1)
discrete_intervals = discretise(normal_dist, intervals)
# Each interval gets probability mass according to the normal distribution
```
"""
function discretise(dist::Distributions.UnivariateDistribution, interval::AbstractVector{IntervalArithmetic.Interval{X}}) where X <: Real
    
    xs_values = vcat(IntervalArithmetic.inf.(interval), [IntervalArithmetic.sup(interval[end])])

    probability = diff(pseudo_cdf.(dist, xs_values))

    #set to sum to one
    probability /= sum(probability)

    return Distributions.DiscreteNonParametric(interval, probability; check_args=false)
end