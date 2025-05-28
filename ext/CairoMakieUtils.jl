"""
   CairoMakie.lines!(ax, z::InvariantCircle; N::Integer=100, color=nothing,
                     i_circle::Integer=0, linewidth=1)

Plot the invariant circle `z` on the CairoMakie axis `ax`.

Arguments:
- `ax`: CairoMakie Axis object
- `z`: The circle in R² to be plotted
- `N`: Number of points to plot
- `i_circle`: Which invariant circle of an island to plot. If 0, plot all
- `color`, `linewidth`: see `CairoMakie.lines!`
"""
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

Useful for plotting with invariant circles on the torus. I.e., if F : T×R→T×R,
and one finds an invariant circle of z(θ+τ) = (h∘F)(z(θ)) where h : T×R→R²,
this plots h⁻¹∘z, the original invariant circle.

Arguments:
- `ax`: CairoMakie Axis object
- `z`: The circle in R² to be plotted
- `hinv`: The map to the torus by the real numbers h⁻¹ : R²→T×R
- `N`: Number of points to plot
- `i_circle`: Which invariant circle of an island to plot. If 0, plot all
- `color`, `linewidth`, `label`: see `CairoMakie.lines!`
"""
function lines_periodic!(ax, z::InvariantCircle, hinv::Function;
                         N::Integer=100, color=nothing, i_circle::Integer=0,
                         linewidth=1, label=nothing)
   #
   if i_circle == 0
      for i_p = 1:get_p(z)
          lines_periodic!(ax, z, hinv; N, color, i_circle=i_p, linewidth, label)
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
             CairoMakie.lines!(ax, hinvs[1,ind], hinvs[2,ind]; linewidth, label)
          else
             CairoMakie.lines!(ax, hinvs[1,ind], hinvs[2,ind]; linewidth, color, label)
          end
      end
   end
end



## Kernel plots
"""
    plot_on_grid(x::AbstractVector, y::AbstractVector, k::KernelLabel;
                 kwargs...)

Create a filled contour plot of the kernel label `k` on the `x` × `y` grid.

kwargs:
- `balance=true`: If true, makes the maximum value of the color scale equal to
  the minimum
- `size`, `fontsize`: See `CairoMakie.Figure`
- `xlabel`, `ylabel`, `title`: See `CairoMakie.Axis`
- `levels`, `linewidth`: See `CairoMakie.contour!`
- `clabel`: See `CairoMakie.Colorbar`

Output:
- `f`: The CairoMakie Figure
- `f_grid`: The data used to make the figure
"""
function plot_on_grid(x::AbstractVector, y::AbstractVector, k::KernelLabel;
                      xlabel="", ylabel="", title="", levels=5, linewidth=2,
                      size = (800, 800), clabel="", fontsize=25,
                      balance=true)
   Nx = length(x);
   Ny = length(y);

   grid = zeros(2, Nx, Ny);
   for ii = 1:Nx, jj = 1:Ny
      grid[:, ii, jj] = [x[ii], y[jj]];
   end
   grid = reshape(grid, 2, Nx*Ny);
   f_grid = reshape(evaluate(k, grid), Nx, Ny);
   f_max = maximum(abs.(f_grid));

   colorrange = (-f_max, f_max)

   if !balance
      colorrange = (minimum(f_grid), maximum(f_grid));
   end

   f = CairoMakie.Figure(size = size, fontsize=fontsize);
   ax = CairoMakie.Axis(f[1,1], xlabel=xlabel, ylabel=ylabel, title=title)

   colormap = colorschemes[:diverging_linear_bjr_30_55_c53_n256]

   levels = LinRange(colorrange[1], colorrange[2], levels)
   p = CairoMakie.contourf!(x, y, f_grid; levels, colormap);
   CairoMakie.contour!(x, y, f_grid; levels, color=:black, linewidth);
   CairoMakie.xlims!(minimum(x), maximum(x))
   CairoMakie.ylims!(minimum(y), maximum(y))
   CairoMakie.Colorbar(f[1, 2], p, label = clabel)

   return f, f_grid
end


"""
    poincare_plot(xb::AbstractVector, yb::AbstractVector, F::Function,
                  Ninit::Integer, Niter::Integer; size=(800, 800),
                  fontsize=25, xlabel="x", ylabel="y", xlims = nothing,
                  ylims=nothing, markersize=3, title="Poincare Plot")

Create a Poincare plot of a 2D map `F` in a rectangular region `xb`×`yb`.
The Poincare plot uses `Ninit` trajectories of length `Niter`.

Output:
- `f`: The figure
- `xs`: The trajectories used to make the figure
"""
function poincare_plot(xb::AbstractVector, yb::AbstractVector,
                       F::Function, Ninit::Integer, Niter::Integer;
                       size=(800, 800), fontsize=25, xlabel="x",
                       ylabel="y", xlims = nothing, ylims=nothing, markersize=3,
                       title="Poincare Plot")
   s = SobolSeq([xb[1], yb[1]], [xb[2], yb[2]])

   f = CairoMakie.Figure(;size, fontsize)
   ax = CairoMakie.Axis(f[1,1]; xlabel, ylabel, title)
   if xlims == nothing; CairoMakie.xlims!(xb...); else; CairoMakie.xlims!(xlims...); end
   if ylims == nothing; CairoMakie.ylims!(yb...); else; CairoMakie.ylims!(ylims...); end

   colors = distinguishable_colors(Ninit, [RGB(1,1,1), RGB(0,0,0)])

   xs = zeros(2, Niter, Ninit)
   for jj = 1:Ninit
      # Get trajectory
      xs[:, 1, jj] = next!(s)
      for ii = 2:Niter
         xs[:, ii, jj] = F(xs[:, ii-1, jj])
      end

      CairoMakie.scatter!(xs[1,:,jj], xs[2,:,jj]; markersize, color=colors[jj])
   end

   f, xs
end
