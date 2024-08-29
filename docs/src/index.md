# SymplecticMapTools.jl
This is a package devoted to different tools for analyzing (primarily 2D)
symplectic maps. The overall goal of this package is to provide a set of tools
that can be used to robustly characterize orbits and find invariant structures
while being efficient in the number of evaluations of the map. If you have any
algorithms that you would like to see added, pull requests are welcome!

The implemented algorithms in this package include
1. The Birkhoff extrapolation method for finding invariant circles [1]
2. The reproducing kernel Hilbert space based invariant level-set finding method [2]
3. The implementations of the parameterization method for finding invariant
   circles and connecting orbits (see, e.g., [3] for an introduction to the
   parameterization method)
4. Newton's method with line search and BFGS for finding periodic orbits
Examples of how to use the code are found in `examples/`. These examples are
created using `Literate.jl` out of files found in `docs/`, and can be recreated
locally by running `docs/literate_examples.jl`. Additionally, these notebooks
are outputted to markdown and included in the documentation.

[1] M. Ruth and D. Bindel, Finding Birkhoff Averages via Adaptive Filtering, Mar. 2024. arXiv:2403.19003 [math.DS].\
[2]  M. Ruth and D. Bindel, Level Set Learning for Poincaré Plots of Symplectic Maps, Dec. 2023. arXiv:2312.00967 [physics].\
[3] A. Haro, M. Canadell, J.-L. Figueras, A. Luque, and J. M. Mondelo,
The Parameterization Method for Invariant Manifolds: From Rigorous Results to
Effective Computations, vol. 195 of Applied Mathematical Sciences,
Springer International Publishing, Cham, 2016.

## Installation
For installation, add the package using the command
`]add SymplecticMapTools` in the Julia REPL.
This package uses Requires for plotting.
You must load CairoMakie or Plots separately to use the plotting functionality.

## Development
This is a new package, and things are changing rapidly! Currently, features
involving the InvariantCircle struct, the extrapolation algorithm, and
periodic orbits (1, 3, and 4 above) are mostly tested, with the features related
to (3) still in need of testing. More examples are being added as well.

## Examples
```@contents
Pages = [
   "examples/birkhoff_averaging/birkhoff_averaging.md"
   "examples/extrapolation/extrapolation.md"
   "examples/kernel/kernel.md"
]
```

## Documentation
```@contents
Pages = [
  "lib/Documentation.md"
  "lib/Internal.md"
]
```

## Contact
For any questions, email [mer335@cornell.edu](mailto:mer335@cornell.edu)
