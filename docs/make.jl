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
    modules = [SymplecticMapTools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/maxeruth/SymplecticMapTools.jl.git"
)
