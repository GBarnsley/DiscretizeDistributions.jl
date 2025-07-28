@doc raw"""
# DiscretiseDistributions.jl

A Julia package for discretising continuous and discrete probability distributions into discrete
representations with specified interval structures.

## Overview

This package provides functionality to convert probability distributions into discrete forms by:
- Discretising continuous distributions using fixed intervals or custom boundaries
- Aggregating discrete distributions into coarser interval representations
- Adjusting discrete distribution alignments (centering and right-alignment), by default discretised distributions are left-aligned, mean the value x gets the probability mass of the interval [x, x + interval)

## Exported Functions

### Main Functions
- [`discretise`](@ref): Convert distributions into discrete representations
- [`center_distribution`](@ref): Center discrete distributions within intervals
- [`right_align_distribution`](@ref): Right-align discrete distributions

## Examples

```julia
using Distributions, DiscretiseDistributions

# Discretise a continuous distribution
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

## Mathematical Details

### Continuous Distributions
For continuous distributions, discretisation computes probability masses using the cumulative
distribution function (CDF):
```math
P(X ∈ [a_i, a_{i+1})) = F(a_{i+1}) - F(a_i)
```

### Discrete Distributions  
For discrete distributions, probability masses are aggregated over intervals using the
probability mass function (PMF):
```math
P(X ∈ [a_i, a_{i+1})) = ∑_{k=⌊a_i⌋}^{⌊a_{i+1}⌋-1} P(X = k)
```

All resulting discrete distributions are normalized to ensure probabilities sum to 1.
"""
module DiscretiseDistributions
    export discretise
    export center_distribution, right_align_distribution
    
    import Distributions, IntervalArithmetic
    import Base: allunique
    
    include("utility.jl")
    include("continous_univariate.jl")
    include("discrete_univariate.jl")
    
    dist = Distributions.Normal(0, 1)
    interval = 0.1
end