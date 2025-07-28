# DiscretiseDistributions.jl

[![Test Status](https://github.com/GBarnsley/DiscretiseDistributions.jl/workflows/CI/badge.svg)](https://github.com/GBarnsley/DiscretiseDistributions.jl/actions)

A Julia package for discretising continuous and discrete probability distributions into discrete representations with specified interval structures.

## Overview

DiscretiseDistributions.jl provides functionality to convert probability distributions into discrete forms by:

- **Discretising continuous distributions** using fixed intervals or custom boundaries
- **Aggregating discrete distributions** into coarser interval representations  
- **Adjusting discrete distribution alignments** (centering and right-alignment)

## Installation

```julia
using Pkg
Pkg.add("DiscretiseDistributions")
```

## Quick Start

```julia
using Distributions, DiscretiseDistributions

# Discretise a continuous distribution with fixed intervals
normal_dist = Normal(0, 1)
discrete_normal = discretise(normal_dist, 0.1)  # Fixed interval width

# Use custom intervals
custom_intervals = [-3.0, -1.0, 0.0, 1.0, 3.0]
discrete_custom = discretise(normal_dist, custom_intervals)

# Discretise a discrete distribution
poisson_dist = Poisson(3.0)
discrete_poisson = discretise(poisson_dist, 2)  # Group into intervals of width 2

# Adjust distribution alignment
centered = center_distribution(discrete_normal)
right_aligned = right_align_distribution(discrete_normal, 0.5)
```

## Exported Functions

### Main Functions

#### `discretise(dist, interval; kwargs...)`

Convert distributions into discrete representations.

**For Continuous Distributions:**
- `discretise(dist::ContinuousUnivariateDistribution, interval::Real; min_quantile=0.001, max_quantile=0.999)`
- `discretise(dist::ContinuousUnivariateDistribution, interval::AbstractVector)`

**For Discrete Distributions:**
- `discretise(dist::DiscreteUnivariateDistribution, interval::Real; min_quantile=0.001, max_quantile=0.999)`  
- `discretise(dist::DiscreteUnivariateDistribution, interval::AbstractVector)`

**Examples:**
```julia
# Fixed interval discretisation
exponential_dist = Exponential(2.0)
discrete_exp = discretise(exponential_dist, 0.5)

# Custom intervals
custom_intervals = [0.0, 1.0, 2.0, 5.0, 10.0]
discrete_custom = discretise(exponential_dist, custom_intervals)

# Control bounds with quantiles
discrete_bounded = discretise(exponential_dist, 0.1; min_quantile=0.01, max_quantile=0.95)
```

#### `center_distribution(dist, interval=nothing)`

Center discrete distributions within intervals.

- `center_distribution(dist::DiscreteNonParametric)` - Center between consecutive support points
- `center_distribution(dist::DiscreteNonParametric, interval)` - Shift all points by `interval/2`

**Examples:**
```julia
dist = DiscreteNonParametric([1.0, 2.0, 3.0, 4.0], [0.25, 0.25, 0.25, 0.25])

# Center between intervals
centered = center_distribution(dist)  # Support: [1.5, 2.5, 3.5]

# Center by specific amount  
centered_shifted = center_distribution(dist, 1.0)  # Support: [1.5, 2.5, 3.5, 4.5]
```

#### `right_align_distribution(dist, interval=nothing)`

Right-align discrete distributions.

- `right_align_distribution(dist::DiscreteNonParametric)` - Shift probabilities to next support point
- `right_align_distribution(dist::DiscreteNonParametric, interval)` - Shift all points by `interval`

**Examples:**
```julia
dist = DiscreteNonParametric([1.0, 2.0, 3.0, 4.0], [0.2, 0.3, 0.3, 0.2])

# Right-align to next support points
right_aligned = right_align_distribution(dist)  # Support: [2.0, 3.0, 4.0]

# Shift by specific amount
right_shifted = right_align_distribution(dist, 0.5)  # Support: [1.5, 2.5, 3.5, 4.5]
```

## Mathematical Details

### Continuous Distributions

For continuous distributions, discretisation computes probability masses using the cumulative distribution function (CDF):

```math
P(X ∈ [a_i, a_{i+1})) = F(a_{i+1}) - F(a_i)
```

where `F(x)` is the CDF of the continuous distribution.

### Discrete Distributions  

For discrete distributions, probability masses are aggregated over intervals using the probability mass function (PMF):

```math
P(X ∈ [a_i, a_{i+1})) = ∑_{k=⌊a_i⌋}^{⌊a_{i+1}⌋-1} P(X = k)
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

## Performance Considerations

- **Pre-allocate intervals**: For repeated discretisations, pre-define interval vectors
- **Choose appropriate quantiles**: For unbounded distributions, adjust quantiles based on your accuracy requirements
- **Interval resolution**: Balance between accuracy and computational efficiency

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Related Packages

- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) - The foundation for probability distributions in Julia
- [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl) - Basic statistics and probability functions

## Citation

If you use DiscretiseDistributions.jl in your research, please consider citing:

```bibtex
@software{DiscretiseDistributions.jl,
  author = {Greg Barnsley},
  title = {DiscretiseDistributions.jl: A Julia package for discretising probability distributions},
  url = {https://github.com/GBarnsley/DiscretiseDistributions.jl},
  year = {2025}
}
```
