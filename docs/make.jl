using Documenter, DiscretizeDistributions

makedocs(
    sitename = "Discretize Distributions",
    pages = [
        "index.md",
        "double_censoring.md"
    ]
)

deploydocs(
    repo = "github.com/GBarnsley/DiscretizeDistributions.jl.git",
)
