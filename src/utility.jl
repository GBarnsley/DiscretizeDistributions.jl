function center_distribution(dist::Distributions.DiscreteNonParametric)
    xs = dist.support[1:(end - 1)] .+ (diff(dist.support) ./ 2)
    return DiscreteNonParametric(xs, dist.p[1:(end - 1)])
end

function center_distribution(dist::Distributions.DiscreteNonParametric, interval)
    xs = dist.support .+ (interval / 2)
    return DiscreteNonParametric(xs, dist.p)
end

function right_align_distribution(dist::Distributions.DiscreteNonParametric)
    return DiscreteNonParametric(dist.support[2:end], dist.p[1:(end - 1)])
end

function right_align_distribution(dist::Distributions.DiscreteNonParametric, interval)
    xs = dist.support .+ (interval)
    return DiscreteNonParametric(xs, dist.p)
end