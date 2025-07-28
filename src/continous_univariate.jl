function define_intervals(dist::Distributions.ContinuousUnivariateDistribution, interval, min_quantile, max_quantile)
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

    return IntervalArithmetic.interval(lower_bound, upper_bound)
end

#a bit of piracy
function allunique(A::Vector{IntervalArithmetic.Interval{X}}) where X <: Real
    if length(A) < 32
        _indexed_allunique(A)
    elseif Base.OrderStyle(eltype(A)) === Base.Ordered()
        a1, rest1 = Iterators.peel(A)::Tuple{Any,Any}
        a2, rest = Iterators.peel(rest1)::Tuple{Any,Any}
        if !IntervalArithmetic.isequal_interval(a1, a2)
            compare = IntervalArithmetic.strictprecedes(a1, a2) ?  IntervalArithmetic.strictprecedes : (a,b) ->  IntervalArithmetic.strictprecedes(b,a)
            for a in rest
                if compare(a2, a)
                    a2 = a
                elseif IntervalArithmetic.isequal_interval(a2, a)
                    return false
                else
                    return Base._hashed_allunique(A)
                end
            end
        else # isequal(a1, a2)
            return false
        end
        return true
    else
        Base._hashed_allunique(A)
    end
end

function discretise_univariate_continuous(dist, xs)
    xs_values = vcat(IntervalArithmetic.inf.(xs), [IntervalArithmetic.sup(xs[end])])
    
    probability = diff(map(Base.Fix1(Distributions.cdf, dist), xs_values))

    #set to sum to one
    probability /= sum(probability)
    allunique(xs)
    return Distributions.DiscreteNonParametric(xs, probability)
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

    xs = define_intervals(dist, interval, min_quantile, max_quantile)

    return discretise_univariate_continuous(dist, xs)
end

@doc """
    discretise(dist::Distributions.ContinuousUnivariateDistribution, interval::AbstractVector)

Discretise a continuous univariate distribution using custom interval boundaries.

This function converts a continuous distribution into a discrete one using user-specified
interval boundaries. The probability mass in each interval is computed using the cumulative
distribution function (CDF).

# Arguments
- `dist::Distributions.ContinuousUnivariateDistribution`: The continuous distribution to discretise
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
```
"""
function discretise(dist::Distributions.ContinuousUnivariateDistribution, interval::AbstractVector)

    xs = sort(interval)

    return discretise_univariate_continuous(dist, xs)
end