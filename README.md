# SymplecticMapTools.jl
This is a package devoted to different tools for analyzing (primarily 2D) symplectic maps. The implemented algorithms include
1. The Birkhoff extrapolation method for finding invariant circles [1]
2. The reproducing kernel Hilbert space based invariant level-set finding method [2]
3. The parameterization method for finding invariant circles and connecting orbits [3]


[1] In preparation.\
[2] In preparation.\
[3] A. Haro, M. Canadell, J.-L. Figueras, A. Luque, and J. M. Mondelo, The Parameterization Method for Invariant Manifolds: From Rigorous Results to Effective Computations, vol. 195 of Applied Mathematical Sciences, Springer International Publishing, Cham, 2016.

# Installation
For installation, clone this repository and add the path using `add path/to/SymplecticMapTools.jl` in the Julia REPL.
This package uses Requires for plotting. You must load CairoMakie or Plots separately to use the plotting functionality.

# Contact
For any questions, email [mer335@cornell.edu](mailto:mer335@cornell.edu)

[![Build Status](https://github.com/maxeruth/SymplecticMapTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/maxeruth/SymplecticMapTools.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/gh/maxeruth/SymplecticMapTools.jl/branch/main/graph/badge.svg?token=5L40XV3NZ8)](https://codecov.io/gh/maxeruth/SymplecticMapTools.jl)
