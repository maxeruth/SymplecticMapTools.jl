using Documenter
using SymplecticMapTools
using Pkg; Pkg.add("Plots"); Pkg.add("CairoMakie")
using Plots
using CairoMakie

# include("./literate_examples.jl")

makedocs(
    sitename = "SymplecticMapTools",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [SymplecticMapTools],
    pages = [
        "SymplecticMapTools.jl" => "index.md",
        "Examples" => [
            "Birkhoff Averages" => "examples/birkhoff_averaging/birkhoff_averaging.md",
            "Birkhoff Extrapolation" => "examples/extrapolation/extrapolation.md",
            "Approximately Invariant Kernel Functions" => "examples/kernel/kernel.md"
        ],
        "Documentation" => "lib/Documentation.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/maxeruth/SymplecticMapTools.jl.git"
)
