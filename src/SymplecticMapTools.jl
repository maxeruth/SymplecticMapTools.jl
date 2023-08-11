# TODO: Add Island.jl and InvariantContinuation.jl
# TODO: LoopVectorization can dramatically speed up stuff in KernelLabel.jl
#       This might be worth a profile
# TODO: I wrote a bunch of windowing functions in KernelLabel.jl? Maybe I should
#       include those

module SymplecticMapTools

## Include packages
using QuadGK
using LinearAlgebra
using FFTW
using LinearOperators
using SparseArrays
using Arpack
using IterativeSolvers
using Polynomials

# using LoopVectorization

## Overloaded functions
# LinearAlgebra.norm, Base.setindex!, Base.getindex, Base.similar, Base.length,
# ???Base.(), Base.size

## Files to include
include("./InvariantCircles/InvariantCircles.jl")
include("./InvariantCircles/ConnectingOrbit.jl")
include("./LabelFunction/KernelLabel.jl")
include("./Birkhoff/BirkhoffAveraging.jl")
include("./Examples/Examples.jl")



## Export functions (comments give where they first show up)
# InvariantCircles.jl, FourierCircle.jl
export InvariantCircle, FourierCircle, get_Na, get_p, get_a0, set_a0!, get_Am,
       set_Am!, get_τ, set_τ!, average_radius, evaluate, deval, deriv, area,
       shifted_eval, get_circle_residual, gn_circle
# ConnectingOrbit.jl
export ConnectingOrbit, get_am, set_am!, linear_initial_connecting_orbit,
       gn_connecting!
# kernels.jl, KernelLabel.jl
export KernelLabel, get_matrix, kernel_sample_F, kernel_eigs, kernel_bvp,
       get_energies, kernel_birkhoff
# Examples.jl
export standard_map_F, standard_map_FJ, polar_map
# BirkhoffAveraging.jl, ContinuedFractions.jl, MPE.jl
export vector_mpe_backslash, vector_mpe_iterative, ContFrac, big_cont_frac_eval,
       partial_frac, denoms, big_partial_frac, big_denoms, wba_weight,
       birkhoff_extrapolation, adaptive_birkhoff_extrapolation, sum_stats,
       get_sum_ave, get_circle_info

## Extensions require Julia > 1.9
# if !isdefined(Base, :get_extension)
      include("../ext/SymplecticMapToolsPlotsExt.jl")
      include("../ext/SymplecticMapToolsCairoMakieExt.jl")
# end

end
