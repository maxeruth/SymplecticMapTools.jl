using Documenter
using SymplecticMapTools
using Plots
using CairoMakie

include("./literate_examples.jl")

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
#=deploydocs(
    repo = "<repository url>"
)=#
