# Discretize Distributions.jl

A Julia package for converting continuous and discrete probability distributions into discrete representations with interval-based support using `IntervalArithmetic.jl`.

The package provides functions to discretize univariate distributions into `DiscreteNonParametric` distributions where the support consists of `IntervalArithmetic.Interval` objects. Each interval `[a, b]` represents a probability mass over that range, computed using the cumulative distribution function (CDF) for continuous distributions or aggregated probability mass function (PMF) for discrete distributions.

## Limitations

- **Finite support**: Infinite distributions are truncated using quantile bounds (default 0.1% and 99.9%)
- **Discrete distribution quirks**: Discretizing already-discrete distributions has some limitations and edge cases
- **Non-integer discrete values**: Discrete distributions with non-integer support may behave unexpectedly
- **Numeric means**: Not all distributions have exact numeric means (i.e. truncated Gamma), these are needed for the `:unbiased` method so a backup numeric mean is calculated where possible using the trapezoid rule with `trapezoid_points`.

## Future Work

- Develop better warnings for incompatible distributions
- Support for multivariate distributions

## API Overview

The package provides three main `discretize` methods:

1. **Fixed intervals**: `discretize(dist, interval_width)` - Creates uniform intervals of specified width
2. **Custom boundaries**: `discretize(dist, boundaries)` - Uses custom interval boundaries  
3. **Pre-constructed intervals**: `discretize(dist, intervals)` - Uses pre-built `Interval` objects

All methods return a `DiscreteNonParametric` distribution with support determined by the `method` parameter.

### Method Parameter

The `discretize` functions accept a `method` parameter that controls the output format:

- `:interval` (default): Returns `IntervalArithmetic.Interval` objects as support points
- `:left_aligned`: Returns left endpoints of intervals as point masses
- `:centred`: Returns interval midpoints as point masses
- `:right_aligned`: Returns right endpoints of intervals as point masses
- `:unbiased`: Returns mean-preserving point masses (requires equal interval widths)

```julia
normal_dist = Normal(0, 1)

# Different output methods
intervals = discretize(normal_dist, 0.5; method=:interval)        # Interval objects
left_points = discretize(normal_dist, 0.5; method=:left_aligned)  # Left endpoints
center_points = discretize(normal_dist, 0.5; method=:centred)     # Midpoints
right_points = discretize(normal_dist, 0.5; method=:right_aligned) # Right endpoints
```

#### Unbiased Method

The `:unbiased` method provides mean-preserving discretization designed to minimize the difference between the original distribution's mean and the discretized distribution's mean. This is an implementation from the `discretize` function in the R package `actuar`.

```julia
# Unbiased discretization - preserves mean
normal_dist = Normal(2.0, 1.0)
unbiased_discrete = discretize(normal_dist, 0.2; method=:unbiased)

# Compare means
println("Original mean: ", mean(normal_dist))      # 2.0
println("Unbiased mean: ", mean(unbiased_discrete)) # ≈ 2.0
println("Centered mean: ", mean(discretize(normal_dist, 0.2; method=:centred)))
```
Both preserve the mean but the unbiased gives more control, supporting all values between [min, min + interval, ..., max] or [lower_quantile, lower_quantile + interval, ..., upper_quantile], where as centred by necessity supports [min + interval/2, min + 3*interval/2, ..., max - interval/2].

### Working with Results

```julia
using Distributions, DiscretizeDistributions, IntervalArithmetic

# Discretize a normal distribution
normal_dist = Normal(0, 1)
interval_dist = discretize(normal_dist, 0.5)

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
discrete_normal = discretize(normal_dist, 0.2; min_quantile=0.005, max_quantile=0.995)

# Exponential distribution - unbounded above
exp_dist = Exponential(1.0)  
discrete_exp = discretize(exp_dist, 0.1; max_quantile=0.99)

# Result includes infinite tail intervals
support(discrete_exp)  # [..., interval(4.5, 5.0), interval(5.0, ∞)]
```

### Custom Interval Structures

Create non-uniform discretisations with custom boundaries:

```julia
# Fine resolution near zero, coarser elsewhere
custom_boundaries = [-5.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
discrete_custom = discretize(Normal(0, 1), custom_boundaries)

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

# Discretize using these intervals
normal_dist = Normal(0, 1)
discrete_custom = discretize(normal_dist, intervals)
```

```@docs
discretize
left_align_distribution
centred_distribution
right_align_distribution
```