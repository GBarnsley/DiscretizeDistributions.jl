# Discretise Distributions.jl Distributions

A set of functions for converting continuous and discrete probability distributions into discrete representations with specified interval structures.
In essence this outputs a `DiscreteNonParametric` with support `x` representing `[x, x + interval)` with probability mass `p` where `p` is either the area of the *pdf* over the interval or the sum of the *pmf* over the interval (assuming the distribution is uniform between its discrete values), standardized to sum to 1.

There are a few limitations:
- The support of the discretized distribution is finite, so infinite distributions are truncated to a finite based on set quantiles (default 0.1% and 99.9% quantiles)
- The support for discretizing a discrete function is weird and limited. In general, you should just do this manually, these are just for my personal convenience. Currently, you can attempt this on distributions where it makes no sense, i.e. a categorical or Bernoulli etc. without error. Also, discrete distributions that support non-integer values are also not supported will still be converted.

Future work:
- Integrate `IntervalArithmetic.jl`
- Add support for intervals with infinite support, this will have to wait since `IntervalArithmetic.jl` is purposefully not interoperable with `Distributions.jl`, though `bareintervals` maybe a solution. This will also simplify the API since all information about the intervals would be contained in the output distribution.
- Develop a method of warnings for attempts are discretizing distributions that are not compatible, might not be possible.
- Support for multivariate distributions?

## Mathematical Details

### Continuous Distributions

For continuous distributions, discretisation computes probability masses using the cumulative distribution function (CDF):

```math
P(X' ∈ [a_i, a_{i+1})) = F(a_{i+1}) - F(a_i)
```

where `F(x)` is the CDF of the continuous distribution `X`.

### Discrete Distributions  

For discrete distributions, probability masses are aggregated over intervals using the probability mass function (PMF):

```math
P(X' ∈ [a_i, a_{i+1})) = ∑_{k=⌈a_i⌉}^{⌊a_{i+1}⌋-1} P(X = k) + (P(X = ⌊a_i⌋) × (⌈a_i⌉ - a_i)) + (P(X = ⌊a_{i+1}⌋) × (a_{i+1} - ⌊a_{i+1}⌋))
```

All resulting discrete distributions are normalized to ensure probabilities sum to 1.

## Advanced Usage

### Handling Unbounded Distributions

For distributions with infinite support, use quantile bounds:

```julia
# Normal distribution - unbounded in both directions
normal_dist = Normal(0, 1)
discrete_normal = discretise(normal_dist, 0.2; min_quantile=0.005, max_quantile=0.995)

# Exponential distribution - unbounded above
exp_dist = Exponential(1.0)
discrete_exp = discretise(exp_dist, 0.1; max_quantile=0.99)
```

### Custom Interval Structures

Create non-uniform discretisations:

```julia
# Fine resolution near zero, coarser elsewhere
custom_intervals = [-5.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
discrete_custom = discretise(Normal(0, 1), custom_intervals)
```

### Working with Distribution Alignments

```julia
# Start with a discretised distribution
dist = discretise(Normal(0, 1), 0.5)

# Center the intervals
centered = center_distribution(dist)

# Right-align for cumulative-like interpretation
right_aligned = right_align_distribution(dist)

# Shift by custom amounts
shifted = right_align_distribution(dist, 0.25)
```
## API Reference

```@docs
discretise
right_align_distribution
center_distribution
```