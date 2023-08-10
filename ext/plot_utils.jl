## Invariant Circles
"""
   Plots.plot(z::InvariantCircle; kwargs...)

Creates a plot of the invariant circle `z`. See `Plots.plot!(z::InvariantCircle)`
for a list of keyword arguments.
"""
function Plots.plot(z::InvariantCircle; N::Integer=100, label=nothing,
                    color=nothing, i_circle=0, linewidth=1, linestyle=:solid)
   if i_circle == 0
      p = Plots.plot(z; N, label, color, i_circle=1, linewidth, linestyle);
      for i_P = 2:get_p(z)
         Plots.plot!(p, z; N, label, color, i_circle=i_P, linestyle);
      end
      return p;
   end

   θ = (0:N)*(2π/N);
   zs = z(θ, i_circle);

   if color==nothing
      return Plots.plot(zs[1,:], zs[2,:]; label, linewidth, linestyle);
   end

   return Plots.plot(zs[1,:], zs[2,:]; label, color, linewidth, linestyle)
end

"""
   Plots.plot(p::Plots.Plot, z::InvariantCircle; kwargs...)

Plots the invariant circle z on p.

kwargs:
- `N::Integer`: Number of ``\theta`` points used to plot each circle
- `i_circle::Integer`: Which period of an island chain to plot. Default value
  of `0` plots all circles of an island chain.
- `label`, `color`, `linewidth`, `linestyle`: see Plots.jl
"""
function Plots.plot!(p::Plots.Plot, z::InvariantCircle; N::Integer=100,
                     i_circle::Integer=0,, label=nothing, color=nothing,
                     linewidth=1, linestyle=:solid)
   if i_circle == 0
      for i_p = 1:get_p(z)
         Plots.plot!(p, z; N, label, color, i_circle=i_p, linewidth, linestyle);
      end
   else
      θ = (0:N)*(2π/N);
      zs = z(θ, i_circle)
      if color==nothing
         Plots.plot!(p, zs[1,:], zs[2,:]; label, linewidth, linestyle);
      else
         Plots.plot!(p, zs[1,:], zs[2,:]; label, color, linewidth, linestyle)
      end
   end

   return p;
end

"""
   parametric_plot(z::InvariantCircle; N::Integer=200, i_circle::Integer=1,
                   linewidth=1, linestyle=:solid, label1="x", label2="y",
                   xlabel="θ", plot_min_dθs=true, markersize=5)

Parametric plot of an invariant circle.

kwargs:
- `N::Integer`: Number of ``\theta`` points used to plot each circle
- `i_circle::Integer`: Which period of an island chain to plot
"""
function parametric_plot(z::InvariantCircle; N::Integer=200, i_circle::Integer=1,
                         linewidth=1, linestyle=:solid, label1="x", label2="y",
                         xlabel="θ", plot_min_dθs=true, markersize=5)
   #
   θ = (0:N)*(2π/N);
   zs = z(θ, i_circle);
   w = 5

   p = Plots.plot(θ, zs[1,:]; label=label1, color=:blue, linestyle, linewidth)
   Plots.plot!(θ, zs[2,:]; label=label2, color=:red, linestyle, linewidth)
   Plots.xlabel!(xlabel);

   if !plot_min_dθs
      return p
   end

   dz = [norm(deval(z, θi)) for θi in θ]
   inds = argminima(dz, w);

   Plots.scatter!(θ[inds], zs[1, inds]; label=false, color=:black, markersize)
   Plots.scatter!(θ[inds], zs[2, inds]; label=false, color=:black, markersize)

   return p, θ[inds], zs[:, inds]
end

## Connecting Orbits
"""
   Plots.plot(c::ConnectingOrbit; kwargs...)

Creates a plot of the connecting orbit `c`. See `Plots.plot!(c::ConnectingOrbit)`
for a list of keyword arguments.
"""
function Plots.plot(c::ConnectingOrbit; N::Integer=50, i_circle=0,
                    label=nothing, color=nothing, linewidth=1)
    p = Plots.plot();
    if i_circle == 0
        for i_P = 1:get_p(c)
            Plots.plot!(p, c; N, label, color, i_circle=i_P, linewidth);
        end
        return p;
    else
        Plots.plot!(p, c; N, label, color, i_circle, linewidth);
    end

   return p
end

"""
   Plots.plot!(p::Plots.Plot, c::ConnectingOrbit; kwargs...)

Creates a plot of the connecting orbit `c`.

kwargs:
- `N::Integer`: Number of points used to plot each orbit
- `i_circle::Integer`: Which period of a connecting orbit chain to plot.
  Default value of `0` plots all connecting orbits of a island chain.
- `label`, `color`, `linewidth`: see Plots.jl
"""
function Plots.plot!(p::Plots.Plot, c::ConnectingOrbit; N::Integer=50, label=nothing,
                     color=nothing, i_circle::Integer=0, linewidth=1)
   if i_circle == 0
      for i_p = 1:get_p(c)
         Plots.plot!(p, c; N, label, color, i_circle=i_p, linewidth);
      end
   else
      s = LinRange(-1.0, 1.0, N);
      zs = c(s, i_circle)
      if color==nothing
         Plots.plot!(p, zs[1,:], zs[2,:]; label, linewidth);
      else
         Plots.plot!(p, zs[1,:], zs[2,:]; label, color, linewidth)
      end
   end

   return p;
end
