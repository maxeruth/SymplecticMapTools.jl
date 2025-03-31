# TODO: Add Island.jl and InvariantContinuation.jl
# TODO: LoopVectorization can dramatically speed up stuff in KernelLabel.jl
#       This might be worth a profile
# TODO: I wrote a bunch of windowing functions in KernelLabel.jl? Maybe I should
#       include those
# TODO: Tests
# TODO: Error handling

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
using Requires
using Sobol
using Optim
using ColorSchemes
using Colors
using LLLplus
using Distributions
using JLD2

# using LoopVectorization

## Overloaded functions
# LinearAlgebra.norm, Base.setindex!, Base.getindex, Base.similar, Base.length,
# ???Base.(), Base.size

## Files to include
include("./PeriodicOrbits/PeriodicOrbits.jl")
include("./InvariantCircles/InvariantCircles.jl")
include("./InvariantTori/FourierTorus.jl")
include("./InvariantCircles/ConnectingOrbit.jl")
include("./LabelFunction/KernelLabel.jl")
include("./Birkhoff/BirkhoffAveraging.jl")
include("./LyapunovExponents/LyapunovExponents.jl")
include("./Examples/Examples.jl")


function __init__()
      @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
            include("../ext/PlotsUtils.jl")
            export parametric_plot
      end

      @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" begin
            include("../ext/CairoMakieUtils.jl")
            export lines_periodic!, plot_on_grid, poincare_plot
      end
end

## Export functions (comments give where they first show up)
# PeriodicOrbits.jl
export BFGS_periodic, newton_periodic

# InvariantCircles.jl, FourierCircle.jl
export InvariantCircle, FourierCircle, get_Na, get_p, get_a0, set_a0!, get_Am, set_Am!, get_τ, 
       set_τ!, average_radius, evaluate, deval, deriv, area, circle_linear!, shifted_eval, 
       get_circle_residual, gn_circle

# FourierTorus.jl
export InvariantTorus, FourierTorus, evaluate_on_grid, kam_residual, kam_residual_norm

# ConnectingOrbit.jl
export ConnectingOrbit, get_am, set_am!, linear_initial_connecting_orbit, gn_connecting!

# kernels.jl, KernelLabel.jl
export KernelLabel, get_matrix, kernel_sample_F, window_weight, rectangular_window_weight, 
       kernel_eigs, kernel_bvp, get_energies, kernel_birkhoff

# Examples.jl
export standard_map_F, standard_map_FJ, polar_map

# BirkhoffAveraging.jl, ContinuedFractions.jl, MPE.jl
export vector_mpe_backslash, vector_mpe_iterative, vector_rre_backslash, vector_rre_iterative,
       ContFrac, big_cont_frac_eval, partial_frac, denoms, big_partial_frac, big_denoms, wba_weight,
       weighted_birkhoff_average, BRREsolution, save_rre, load_rre, doubling_birkhoff_average, 
       birkhoff_extrapolation, adaptive_birkhoff_extrapolation, get_w0, get_w0!, get_torus, 
       sum_stats, get_sum_ave, get_circle_info, adaptive_get_torus, adaptive_get_torus!

# LyapunovExponents.jl
export lyapunov_exponent, lyapunov_exponents


# using CairoMakie
# include("../ext/CairoMakieUtils.jl")
# export lines_periodic!, plot_on_grid, poincare_plot

## Extensions require Julia > 1.9
# if !isdefined(Base, :get_extension)

# end

end
