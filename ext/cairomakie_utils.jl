## Invariant Circles
function CairoMakie.lines!(ax, z::InvariantCircle; N::Integer=100,
                           color=nothing, i_circle::Integer=0,
                           linewidth=1)
   #
   if i_circle == 0
      for i_p = 1:get_p(z)
         CairoMakie.lines!(ax, z; N, color, i_circle=i_p, linewidth)
      end
   else
      θ = (0:N)*(2π/N);
      zs = z(θ, i_circle);
      if color==nothing
         CairoMakie.lines!(ax, zs[1,:], zs[2,:]; linewidth)
      else
         CairoMakie.lines!(ax, zs[1,:], zs[2,:]; linewidth, color)
      end
   end
end

"""
   lines_periodic!(ax, z::InvariantCircle, hinv::Function; N::Integer=100,
                   color=nothing, i_circle::Integer=0, linewidth=1)

Useful for plotting with invariant circles on the torus. I.e., if F : 𝕋×ℝ→𝕋×ℝ,
and one finds an invariant circle of z(θ+τ) = (h∘F)(z(θ)) where h : 𝕋×ℝ→ℝ²,
this plots h⁻¹∘z, the original invariant circle.

Arguments:
- `ax`: CairoMakie Axis object
- `z`: The circle in ℝ² to be plotted
- `hinv`: The map to the torus by the real numbers h⁻¹ : ℝ²→𝕋×ℝ
- `N`: Number of points to plot
- `i_circle`: Which invariant circle of an island to plot. If 0, plot all
- `color`, `linewidth`: see `CairoMakie.lines!`
"""
function lines_periodic!(ax, z::InvariantCircle, hinv::Function;
                         N::Integer=100, color=nothing, i_circle::Integer=0,
                         linewidth=1)
   #
   if i_circle == 0
      for i_p = 1:get_p(z)
          lines_periodic!(ax, z, hinv; N, color, i_circle=i_p, linewidth)
      end
   else
      θ = (0:N)*(2π/N);
      zs = z(θ, i_circle);

      hinvs = similar(zs);
      for (ii, z) in enumerate(eachcol(zs))
          hinvs[:, ii] = hinv(z);
      end
      diffnorms = [norm(diffcol) for diffcol in eachcol(diff(hinvs, dims=2))]
      h_ave = (hinvs*ones(N+1)) ./ (N+1)
      hscale = maximum([norm(hcol-h_ave) for hcol in eachcol(hinvs)])


      jumps = vcat(0, (1:N)[diffnorms .> 0.1*hscale], N+1)
      # jumps = [0, N+1]
      for ii = 1:length(jumps)-1
          ind = jumps[ii]+1:jumps[ii+1]
          if color==nothing
             CairoMakie.lines!(ax, hinvs[1,ind], hinvs[2,ind]; linewidth)
          else
             CairoMakie.lines!(ax, hinvs[1,ind], hinvs[2,ind]; linewidth, color)
          end
      end
   end
end



## Kernel plots
"""
   plot_on_grid(x::AbstractVector, y::AbstractVector, k::KernelLabel; kwargs...)

Create a contour plot of the kernel label `k` on the `x` × `y` grid.

kwargs:
- `balance=true`: If true, makes the maximum value of the color scale equal to
  the minimum
- `resolution`, `fontsize`: See `CairoMakie.Figure`
- `xlabel`, `ylabel`, `title`: See `CairoMakie.Axis`
- `levels`, `linewidth`: See `CairoMakie.contour!`
- `clabel`: See `CairoMakie.Colorbar`
"""
function plot_on_grid(x::AbstractVector, y::AbstractVector, k::KernelLabel;
                      xlabel="", ylabel="", title="", levels=5, linewidth=2,
                      resolution = (800, 800), clabel="", fontsize=30, balance=true)
    Nx = length(x);
    Ny = length(y);

    grid = zeros(2, Nx, Ny);
    for ii = 1:Nx, jj = 1:Ny
       grid[:, ii, jj] = [x[ii], y[jj]];
    end
    grid = reshape(grid, 2, Nx*Ny);
    f_grid = reshape(eval(k, grid), Nx, Ny);
    f_max = maximum(abs.(f_grid));
    colorrange = (-f_max, f_max)

    if !balance
        colorrange = (minimum(f_grid), maximum(f_grid));
    end

    f = Figure(resolution = resolution, fontsize=fontsize);
    ax = Axis(f[1,1], xlabel=xlabel, ylabel=ylabel, title=title)

    colormap = "RdBu"
    # colormap = range(HSL(colorant"red"), stop=HSL(colorant"blue"), length=15)

    p = CairoMakie.contour!(x, y, f_grid; levels=levels, linewidth=linewidth,
                                        colormap=colormap, colorrange=colorrange);
    Colorbar(f[1, 2], p, label = clabel)

    return f, ax, f_grid
end
