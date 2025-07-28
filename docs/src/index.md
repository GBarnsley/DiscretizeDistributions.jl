# Discretise Distributions.jl

A Julia package for converting continuous and discrete probability distributions into discrete representations with interval-based support using `IntervalArithmetic.jl`.

The package provides functions to discretise univariate distributions into `DiscreteNonParametric` distributions where the support consists of `IntervalArithmetic.Interval` objects. Each interval `[a, b]` represents a probability mass over that range, computed using the cumulative distribution function (CDF) for continuous distributions or aggregated probability mass function (PMF) for discrete distributions.

## Limitations

- **Finite support**: Infinite distributions are truncated using quantile bounds (default 0.1% and 99.9%)
- **Discrete distribution quirks**: Discretizing already-discrete distributions has some limitations and edge cases
- **Non-integer discrete values**: Discrete distributions with non-integer support may behave unexpectedly

## Future Work

- Develop better warnings for incompatible distributions
- Support for multivariate distributions

## API Overview

The package provides three main `discretise` methods:

1. **Fixed intervals**: `discretise(dist, interval_width)` - Creates uniform intervals of specified width
2. **Custom boundaries**: `discretise(dist, boundaries)` - Uses custom interval boundaries  
3. **Pre-constructed intervals**: `discretise(dist, intervals)` - Uses pre-built `Interval` objects

All methods return a `DiscreteNonParametric` distribution with `IntervalArithmetic.Interval` support.

### Working with Results

```julia
using Distributions, DiscretiseDistributions, IntervalArithmetic

# Discretise a normal distribution
normal_dist = Normal(0, 1)
interval_dist = discretise(normal_dist, 0.5)

# The result has interval support
support(interval_dist)  # Vector of Interval{Float64} objects
probs(interval_dist)    # Corresponding probabilities

# Convert to point-based distributions
left_aligned = left_align_distribution(interval_dist)     # Use left endpoints
centered = centred_distribution(interval_dist)            # Use midpoints  
right_aligned = right_align_distribution(interval_dist)   # Use right endpoints
```

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

For distributions with infinite support, control truncation with quantile bounds:

```julia
# Normal distribution - unbounded in both directions  
normal_dist = Normal(0, 1)
discrete_normal = discretise(normal_dist, 0.2; min_quantile=0.005, max_quantile=0.995)

# Exponential distribution - unbounded above
exp_dist = Exponential(1.0)  
discrete_exp = discretise(exp_dist, 0.1; max_quantile=0.99)

# Result includes infinite tail intervals
support(discrete_exp)  # [..., interval(4.5, 5.0), interval(5.0, ∞)]
```

### Custom Interval Structures

Create non-uniform discretisations with custom boundaries:

```julia
# Fine resolution near zero, coarser elsewhere
custom_boundaries = [-5.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
discrete_custom = discretise(Normal(0, 1), custom_boundaries)

# Results in intervals: [(-∞,-5], [-5,-2], [-2,-1], ..., [5,∞)]
length(support(discrete_custom))  # 10 intervals (8 from boundaries + 2 infinite tails)
```

### Working with Pre-constructed Intervals

For advanced use cases, you can provide pre-constructed `IntervalArithmetic.Interval` objects:

```julia
using IntervalArithmetic

# Create custom intervals with specific properties
intervals = [
    interval(-2.0, -1.0),    # Standard interval
    interval(-1.0, 0.0),     # Adjacent interval
    interval(0.0, 2.0),      # Wider interval
    interval(2.0, Inf)       # Semi-infinite interval
]

# Discretise using these intervals
normal_dist = Normal(0, 1)
discrete_custom = discretise(normal_dist, intervals)
```

```@docs
discretise
left_align_distribution
centred_distribution
right_align_distribution
```