using Documenter, DiscretiseDistributions

makedocs(
    sitename="Discretise Distributions.jl Distributions",
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo = "github.com/GBarnsley/DiscretiseDistributions.jl.git",
)