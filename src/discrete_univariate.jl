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
    probability = Distributions.pdf(dist, xs_daily)

    #set to sum to one
    probability /= sum(probability)

    return Distributions.DiscreteNonParametric(xs[probability .> 0], probability[probability .> 0])
end

function discretise(dist::Distributions.DiscreteUnivariateDistribution, interval::Real; min_quantile = 0.001, max_quantile=0.999)

    range = define_bounds(dist, interval, min_quantile, max_quantile)
    
    xs = collect(range * interval)

    return discretise_univariate_discrete(dist, xs)
end

function discretise(dist::Distributions.DiscreteUnivariateDistribution, interval::AbstractVector)

    xs = sort(interval)

    return discretise_univariate_discrete(dist, xs)
end