# SymplecticMapTools.jl

Documentation for SymplecticMapTools.jl

## Contents

```@contents
Pages = ["index.md"]
```

## Index

```@index
Pages = ["index.md"]
```

## Invariant Circles
```@docs
InvariantCircle
evaluate(::InvariantCircle, ::Number)
deval(::InvariantCircle, ::Number)
shifted_eval
get_circle_residual
gn_circle
```

## Fourier Circles
```@docs
FourierCircle
get_Na(::FourierCircle)
get_p(::FourierCircle)
Base.length(::FourierCircle)
get_a0
set_a0!
get_Am
set_Am!
Base.getindex(::FourierCircle, ::Integer, ::Integer)
Base.setindex!(::FourierCircle, ::AbstractArray, ::Integer, ::Integer)
get_τ
set_τ!
average_radius(::FourierCircle)
evaluate(::FourierCircle, ::AbstractVector)
deval(::FourierCircle, ::AbstractVector)
deriv(::FourierCircle)
area
```

## Connecting Orbits
```@docs
ConnectingOrbit
get_Na(::ConnectingOrbit)
get_p(::ConnectingOrbit)
get_am(::ConnectingOrbit, ::Integer)
set_am!(::ConnectingOrbit, ::Integer, ::AbstractArray)
Base.getindex(::ConnectingOrbit, ::Integer, ::Integer)
Base.setindex!(::ConnectingOrbit, ::AbstractArray, ::Integer, ::Integer)
evaluate(::ConnectingOrbit, ::AbstractArray)
evaluate(::ConnectingOrbit, ::Number)
deval(::ConnectingOrbit, ::Number)
area(::ConnectingOrbit)
linear_initial_connecting_orbit
gn_connecting!
```

## Kernel Labels
```@docs
KernelLabel
get_matrix(::KernelLabel, ::AbstractArray)
evaluate(::KernelLabel, ::AbstractArray)
kernel_sample_F
kernel_eigs
kernel_bvp
get_energies
kernel_birkhoff
```

## Continued Fractions
```@docs
ContFrac
big_cont_frac_eval
evaluate(::ContFrac)
partial_frac
denoms
big_partial_frac
big_denoms
```

## Birkhoff Extrapolation
```@docs
vector_mpe_backslash
vector_mpe_iterative
wba_weight
birkhoff_extrapolation
adaptive_birkhoff_extrapolation
sum_stats
get_sum_ave
get_circle_info
```

## Example Maps
```@docs
standard_map_F
standard_map_FJ
polar_map
```

## Plotting Routines: Plots
```@autodocs
Modules = [SymplecticMapTools.SymplecticMapToolsPlotsExt]
```

## Plotting Routines: CairoMakie
```@autodocs
Modules = [SymplecticMapTools.SymplecticMapToolsCairoMakieExt]
```
