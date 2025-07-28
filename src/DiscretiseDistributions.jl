@doc raw"""
# DiscretizeDistributions.jl

A Julia package for discretising continuous and discrete probability distributions into interval-based discrete
representations with flexible alignment options.

## Overview

This package provides functionality to convert probability distributions into discrete forms using intervals by:
- **Discretising distributions using intervals**: Backend creates discrete distributions with `IntervalArithmetic.Interval` support
- **Flexible alignment options**: Convert interval-based distributions to point-based distributions with left, center, or right alignment
- **Fixed or custom intervals**: Use uniform intervals or specify custom interval boundaries

The package follows a two-stage approach:
1. **Backend discretisation**: Creates distributions with interval support using `IntervalArithmetic.jl`
2. **Alignment utilities**: Convert intervals to point values for final output (left-aligned, centered, or right-aligned)

## Exported Functions

### Main Functions
- [`discretize`](@ref): Convert distributions into interval-based discrete representations
- [`left_align_distribution`](@ref): Convert intervals to left-aligned point values
- [`centred_distribution`](@ref): Convert intervals to centered point values  
- [`right_align_distribution`](@ref): Convert intervals to right-aligned point values

## Examples

```julia
using Distributions, DiscretizeDistributions
using IntervalArithmetic

# Discretize a continuous distribution (returns interval-based distribution)
normal_dist = Normal(0, 1)
discrete_intervals = discretize(normal_dist, 0.5)  # Support: intervals like [-1.0, -0.5), [-0.5, 0.0), etc.

# Convert to different alignments
left_aligned = left_align_distribution(discrete_intervals)    # Support: [-1.0, -0.5, 0.0, ...]
centered = centred_distribution(discrete_intervals)          # Support: [-0.75, -0.25, 0.25, ...]  
right_aligned = right_align_distribution(discrete_intervals) # Support: [-0.5, 0.0, 0.5, ...]

# Use custom intervals
custom_intervals = [-3.0, -1.0, 0.0, 1.0, 3.0]
discrete_custom = discretize(normal_dist, custom_intervals)

# Discretize a discrete distribution
poisson_dist = Poisson(3.0)
discrete_poisson = discretize(poisson_dist, 2)  # Group into intervals of width 2
```
"""
module DiscretizeDistributions
    export discretize
    export left_align_distribution, centred_distribution, right_align_distribution
    
    import Distributions, IntervalArithmetic
    
    include("utility.jl")
    include("univariate.jl")
    
end