using Documenter, DiscretizeDistributions

makedocs(
    sitename="Discretize Distributions.jl Distributions",
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo = "github.com/GBarnsley/DiscretizeDistributions.jl.git",
)