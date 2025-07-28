# DiscretiseDistributions.jl

[![Build Status](https://github.com/GBarnsley/DiscretiseDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GBarnsley/DiscretiseDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://GBarnsley.github.io/DiscretiseDistributions.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://GBarnsley.github.io/DiscretiseDistributions.jl/dev)

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

I recommend you use rational numbers to defined your intervals for precision. You can then convert to floating-point numbers if needed.

```julia
using Distributions, DiscretiseDistributions

# Discretise a continuous distribution with fixed intervals
normal_dist = Normal(0.0, 1.0)
discrete_normal = discretise(normal_dist, 1//10)  # Fixed interval width

# Use custom intervals
custom_intervals = [-3.0, -1.0, 0.0, 1.0, 3.0]
discrete_custom = discretise(normal_dist, custom_intervals)

# Due to limitations we have to limit the distribution to finitely bounded intervals
# Currently the bounds for distributions with infinite support are set be the largest/smallest that include the quantiles 0.1% and 99.9% quantiles (customisable)
# You can also set the bound using truncation
maximum(discrete_normal.support) #3
discrete_normal_truncated = discretise(truncated(Normal(), -1.1, 5.13), 1//10)
minimum(discrete_normal_truncated.support) #-1.1
maximum(discrete_normal_truncated.support) #5 (representing 5 - 5.1)
# values between 5.1 and 5.13 are "lost" and do not contribute to probability due to the interval size

# Discretise a discrete distribution
poisson_dist = Poisson(3.0)
discrete_poisson = discretise(poisson_dist, 2)  # Group into intervals of width 2

# Adjust distribution alignment
# convenience functions, you can also do this yourself by accessing .support
centered = center_distribution(discrete_normal)
right_aligned = right_align_distribution(discrete_normal, 0.5)
```

## Related Packages

- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) - The foundation for probability distributions in Julia