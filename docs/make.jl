using Documenter
using SymplecticMapTools

makedocs(
    sitename = "SymplecticMapTools",
    format = Documenter.HTML(),
    modules = [SymplecticMapTools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
