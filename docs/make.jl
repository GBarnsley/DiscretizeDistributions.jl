using Documenter, discretizeDistributions

makedocs(
    sitename="discretize Distributions.jl Distributions",
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo = "github.com/GBarnsley/discretizeDistributions.jl.git",
)