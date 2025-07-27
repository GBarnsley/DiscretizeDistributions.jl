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

function discretise_univariate_continuous(dist, xs)
    probability = diff(Distributions.cdf(dist, xs))

    xs = xs[1:(end-1)]

    #set to sum to one
    probability /= sum(probability)

    return Distributions.DiscreteNonParametric(xs[probability .> 0], probability[probability .> 0])
end

function discretise(dist::Distributions.ContinuousUnivariateDistribution, interval::Real; min_quantile = 0.001, max_quantile=0.999)

    range = define_bounds(dist, interval, min_quantile, max_quantile)

    xs = collect(range * interval)

    return discretise_univariate_continuous(dist, xs)
end

function discretise(dist::Distributions.ContinuousUnivariateDistribution, interval::AbstractVector)

    xs = sort(interval)

    return discretise_univariate_continuous(dist, xs)
end