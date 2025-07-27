module DiscretiseDistributions
    export discretise
    export center_distribution, right_align_distribution
    
    import Distributions
    
    include("utility.jl")
    include("continous_univariate.jl")
    include("discrete_univariate.jl")
    
end