# DiscretizeDistributions.jl

[![Build Status](https://github.com/GBarnsley/DiscretizeDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GBarnsley/DiscretizeDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://GBarnsley.github.io/DiscretizeDistributions.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://GBarnsley.github.io/DiscretizeDistributions.jl/dev)
[![codecov](https://codecov.io/gh/GBarnsley/DiscretizeDistributions.jl/graph/badge.svg)](https://codecov.io/gh/GBarnsley/DiscretizeDistributions.jl)


[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package for discretising continuous and discrete probability distributions into interval-based discrete representations with flexible alignment options.

## Overview

DiscretizeDistributions.jl provides functionality to convert probability distributions into discrete forms using a robust two-stage approach:

1. **Backend discretisation**: Creates discrete distributions with `IntervalArithmetic.Interval` support representing probability masses over intervals
2. **Alignment utilities**: Convert interval-based distributions to point-based distributions for final output, if needed

### Key Features

- **Interval-based backend**: Uses `IntervalArithmetic.jl` for precise interval representation
- **Flexible output formats**: Convert to left-aligned, centered, or right-aligned point distributions
- **Support for both continuous and discrete input distributions**
- **Custom or uniform interval spacing**
- **Automatic handling of infinite bounds using quantiles**

## Installation

```julia
using Pkg
Pkg.add("DiscretizeDistributions")
```

## Quick Start

We recommend using rational numbers for interval definitions to maintain precision:

```julia
using Distributions, IntervalArithmetic
using DiscretizeDistributions

# Discretize a continuous distribution (returns interval-based distribution)
normal_dist = Normal(0.0, 1.0)
discrete_intervals = discretize(normal_dist, 1//10)  # Returns intervals like [0.0, 0.1), [0.1, 0.2), etc.

# Convert to different point-based alignments
left_aligned = left_align_distribution(discrete_intervals)    # Uses interval start points
centered = centred_distribution(discrete_intervals)          # Uses interval midpoints  
right_aligned = right_align_distribution(discrete_intervals) # Uses interval end points

# Use custom intervals
custom_intervals = [-3.0, -1.0, 0.0, 1.0, 3.0]
discrete_custom = discretize(normal_dist, custom_intervals)

# For unbounded distributions, bounds are set using quantiles (customizable)
# Default: 0.1% and 99.9% quantiles
println("Automatic bounds: ", minimum(discrete_intervals.support), " to ", maximum(discrete_intervals.support))

# You can also use truncation for explicit bounds
discrete_truncated = discretize(truncated(Normal(), -1.1, 5.13), 1//10)
println("Truncated bounds: ", minimum(discrete_truncated.support), " to ", maximum(discrete_truncated.support))

# Discretize a discrete distribution
poisson_dist = Poisson(3.0)
discrete_poisson = discretize(poisson_dist, 2)  # Group into intervals of width 2

# The backend maintains interval structure - convert to points as needed
println("Interval support: ", support(discrete_intervals)[1:5])  # Shows first 5 intervals
println("Left-aligned: ", support(left_aligned)[1:5])           # Shows first 5 points
```

## Related Packages

- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) - Foundation for probability distributions in Julia
- [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) - Interval arithmetic backend