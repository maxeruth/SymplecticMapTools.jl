# SymplecticMapTools.jl
This is a package devoted to different tools for analyzing (primarily 2D) symplectic maps. The implemented algorithms include
1. The parameterization method for finding invariant circles [1]
2. The Birkhoff filter diagonalization method for finding invariant circles [2]
3. The reproducing kernel Hilbert space based invariant level-set finding method [3]

[1] A. Haro, M. Canadell, J.-L. Figueras, A. Luque, and J. M. Mondelo, The Parameterization Method for Invariant Manifolds: From Rigorous Results to Effective Computations, vol. 195 of Applied Mathematical Sciences, Springer International Publishing, Cham, 2016.\
[2] In preparation.\
[3] In preparation.

# Installation
For installation, simply use the Julia package manager command `]install SymplecticMapTools`.
This package uses Julia extensions to dynamically load CairoMakie or Plots.
If you do not have >1.9, the package will automatically load the plotting software.

# Contact
For any questions, email [mer335@cornell.edu](mailto:mer335@cornell.edu)

[![Build Status](https://github.com/maxeruth/SymplecticMapTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/maxeruth/SymplecticMapTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
