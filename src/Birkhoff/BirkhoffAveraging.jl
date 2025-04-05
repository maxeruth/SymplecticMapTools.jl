include("./ExtrapolationOperators.jl")
include("./MPE.jl")
include("./ContinuedFractions.jl")

# The distance between two numbers on the torus
function mid_mod(theta::Number)
    mod(theta+0.5, 1.) - 0.5
end

"""
    wba_weight(d, N)

Weights used in the weighted Birkhoff average, returned as a diagonal matrix.
Returns as `Diagonal([w_1, ..., w_1, ..., w_N, ..., w_N])`, where each
coefficient is repeated `d` times (this is used for vector MPE and RRE as a
least-squares weighting)
"""
function wba_weight(d::Integer, N::Integer)
    t = (1:N) ./ (N+1);
    g = exp.(- 1 ./ ( t .* (1 .- t))) .+ 1e-16;
    g = sqrt.(g ./ sum(g))
    return Diagonal(kron(g,ones(d)));
end


"""
    weighted_birkhoff_average(hs::AbstractMatrix)

Finds a weighted Birkhoff average of a sequence of vector observations. The
array input is assumed to be of size d × N, where the average is performed over
the second index.
"""
function weighted_birkhoff_average(hs::AbstractMatrix)
    N = size(hs, 2)
    t = (1:N) ./ (N+1);
    w = exp.(- 1 ./ ( t .* (1 .- t)))
    return (hs*w) ./ sum(w)
end

"""
    weighted_birkhoff_average(hs::AbstractVector)

Finds a weighted Birkhoff average of a sequence of scalar observations.
"""
function weighted_birkhoff_average(hs::AbstractVector)
    weighted_birkhoff_average(hs')
end


"""
    doubling_birkhoff_average(h::Function, F::Function, x0::AbstractVector;
                              tol::Number = 1e-10, T_init::Integer=10,
                              T_max::Integer=320)

Find the weighted Birkhoff ergodic average of an observable function `h` over a
trajectory of the map `F` adaptively.

Arguments:
- `h`: A function from Rᵈ to Rⁿ, where d is the dimension of the state space
  and n is the dimension of the observation
- `F`: The (symplectic) map: a function from Rᵈ to Rᵈ
- `x0`: The initial point of the trajectory in Rᵈ
- `tol`: The tolerance by which convergence is judged. If the average does not
  change by more than `tol` over a doubling, the function returns
- `T_init`: The initial length of the trajectory considered
- `T_max`: The maximum trajectory length considered

Output:
- `ave`: The average of `h` over the trajectory
- `xs`: The trajectory
- `hs`: The value of `h` on the trajectory
- `conv_flag`: `true` iff the averages converged
"""
function doubling_birkhoff_average(h::Function, F::Function, x0::AbstractVector;
                                   tol::Number = 1e-10, T_init::Integer=10,
                                   T_max::Integer=320)
    # Initialize arrays
    T = T_init
    d = length(x0)
    xs = zeros(d, T)
    xs[:, 1] = x0;

    h0 = h(x0)
    d_h = length(h0);
    hs = zeros(d_h, T);
    hs[:, 1] .= h0;

    # Find initial average
    for t = 2:T
        xs[:, t] = F(xs[:, t-1])
        hs[:, t] .= h(xs[:, t])
    end
    ave = weighted_birkhoff_average(hs)

    # Double average size until convergence or errors
    while 2T <= T_max
        xs_prev = xs;
        hs_prev = hs;
        ave_prev = ave

        xs = zeros(d, 2T)
        hs = zeros(d_h, 2T)

        xs[:, 1:T] = xs_prev;
        hs[:, 1:T] = hs_prev;

        for t = T+1:2T
            xs[:, t] = F(xs[:, t-1])
            hs[:, t] .= h(xs[:, t])
        end

        ave = weighted_birkhoff_average(hs)

        if norm(ave-ave_prev) < tol
            return ave, xs, hs, true
        end

        T = 2T
    end


    return ave, xs, hs, false
end

mutable struct BRREsolution
    D::Integer  # Observation Function Dimension
    F::Function # Evaluation Function
    h::Function # Observation Function

    xs::AbstractArray # Trajectory
    hs::AbstractArray # Observations on Trajectory

    K::Integer            # No. Filter Unknowns
    c::AbstractVector     # Filter Coefficients
    h_ave::AbstractVector # Birkhoff Average
    resid_RRE::Number     # RRE Residual

    d::Integer         # Torus Dimension
    w0::AbstractVector # Rotation Vector
    resid_w0::Number   # Rotation Vector Residual

    tor::FourierTorus # Invariant Torus
    resid_tor::Number # Invariant Torus Residual

    function BRREsolution(D::Integer, F::Function, h::Function, xs::AbstractArray, 
                          hs::AbstractArray, K::Integer, c::AbstractVector, h_ave::AbstractVector, 
                          resid_RRE::Number)
        d = 0;
        w0 = Float64[];
        resid_w0 = Inf;

        tor = FourierTorus(D,ones(Integer,d))
        resid_tor = Inf;

        new(D,F,h,xs,hs,K,c,h_ave,resid_RRE,d,w0,resid_w0,tor,resid_tor)
    end
end

"""
    save_rre(file::AbstractString, sol::BRREsolution)

Save an RRE solution object to a file (not including the evaluation and observation
function `F` and `h`). Currently saves using JLD2, but could be 
extended in the future.
"""
function save_rre(file::AbstractString, sol::BRREsolution)
    D = sol.D
    xs = sol.xs
    hs = sol.hs
    K = sol.K
    c = sol.c
    h_ave = sol.h_ave
    resid_RRE = sol.resid_RRE
    d = sol.d
    w0 = sol.w0
    resid_w0 = sol.resid_w0
    tor = sol.tor
    resid_tor = sol.resid_tor

    tor_a = tor.a
    tor_τ = tor.τ
    tor_d, tor_p, tor_Na = tor.sz

    @save file D xs hs K c h_ave resid_RRE d w0 resid_w0 resid_tor tor_a tor_τ tor_d tor_p tor_Na
end
 
"""
    load_rre(file::AbstractString, sol::BRREsolution)

Load an RRE solution object from a file (saved by `save_rre`). 
"""
function load_rre(file::AbstractString)
    @load file D xs hs K c h_ave resid_RRE d w0 resid_w0 resid_tor tor_a tor_τ tor_d tor_p tor_Na
    tor = FourierTorus(tor_d, tor_Na; a=tor_a, p=tor_p, τ=tor_τ)
    sol = BRREsolution(D, (x)->nothing, (x)->nothing, xs, hs, K, c, h_ave, resid_RRE)
    sol.d = d
    sol.w0 = w0
    sol.resid_w0 = resid_w0
    sol.tor = tor
    sol.resid_tor = resid_tor
    
    sol
end

function unorm(hs::AbstractArray)
    us = diff(hs, dims=2)
    us2 = [ui'*ui for ui in eachcol(us)]
    sqrt(weighted_birkhoff_average(us2))
end


"""
    birkhoff_extrapolation(h::Function, F::Function, x0::AbstractVector,
                           N::Integer, K::Integer; iterative::Bool=false,
                           x_prev::Union{AbstractArray,Nothing}=nothing,
                           rre::Bool=false, weighted::Bool=false)

As input, takes an initial point `x0`, a symplectic map `F`, and an observable
`h` (can choose the identity as a default `h = (x)->x`). Then, the method
1. Computes a time series xs[:,n+1] = Fⁿ(x0)
2. Computes observable series hs[:,n] = h(xs[:,n])
3. Performs sequence extrapolation (RRE or MPE) on hs to obtain a model `c` of
    length `2K+1` (with `K` unknowns), the extrapolated value applied to each
    window, and a residual
4. Returns as `c, sums, resid, xs, hs[, history]` where `history` is a
    diagnostic from the iterative solver that only returns when `iterative=true`

Use `x_prev` if you already know part of the sequence, but do not know the whole
thing.
"""
function birkhoff_extrapolation(h::Function, F::Function, x0::AbstractVector,
                                N::Integer, K::Integer; iterative::Bool=false,
                                x_prev::Union{AbstractArray,Nothing}=nothing,
                                rre::Bool=true, ϵ::Number=0.0, weighted::Bool=false,
                                ortho::Bool=true)
    x = deepcopy(x0);
    h0 = h(x);
    D = length(h0);

    @assert N*D ≥ K

    Nx = N+2K+1
    hs = zeros(D, Nx);
    xs = zeros(D, Nx);

    hs[:, 1] = h0;
    xs[:, 1] = x;
    ii_init = 2;
    if x_prev != nothing
        N_prev = size(x_prev, 2)
        xs[:, 1:N_prev] = x_prev;
        for ii = 1:N_prev
            hs[:, ii] = h(xs[:,ii])
        end
        ii_init = N_prev+1
    end
    x = xs[:, ii_init-1]

    for ii = ii_init:Nx
        x = F(x);
        hs[:,ii] = h(x)
        xs[:,ii] = x;
    end
    history = 0

    if rre
        if iterative
            c, sums, resid, history = vector_rre_iterative(hs, K; ϵ, weighted, ortho)
        else
            c, sums, resid = vector_rre_backslash(hs, K; ϵ, weighted, ortho)
        end
    else
        if iterative
            c, sums, resid, history = vector_mpe_iterative(hs, K; ϵ, weighted, ortho)
        else
            c, sums, resid = vector_mpe_backslash(hs, K; ϵ, weighted, ortho)
        end
    end

    h_ave = sums[:,1];

    resid_RRE = norm(resid)/unorm(hs)

    BRREsolution(D, F, h, xs, hs, K, c, h_ave, resid_RRE)
end


"""
    adaptive_birkhoff_extrapolation(h::Function, F::Function,
                    x0::AbstractVector; rtol::Number=1e-10, Kinit = 20,
                    Kmax = 100, Kstride=20, iterative::Bool=false,
                    Nfactor::Integer=5, rre::Bool=true)

Adaptively applies `birkhoff_extrapolation` to find a good enough filter length
`K`, where "good enough" is defined by the `rtol` optional argument.

Arguments:
- `h`: The observable function (often the identity function `h = (x)->x` works)
- `F`: The symplectic map
- `x0`: The initial point of the trajectory
- `rtol`: Required tolerance for convergence (inexact maps often require a
  looser tolerance)
- `Kinit`: The length of the initial filter
- `Kmax`: The maximum allowed filter size
- `Kstride`: The amount `K` increases between applications of
             `birkhoff_extrapolation`
- `iterative`: Whether to use an iterative method to solve the Hankel system
               in the extrapolation step
- `Nfactor`: How rectangular the extrapolation algorithm is. Must be >=1.
- `rre`: Turn to true to use reduced rank extrapolation instead of minimal
         polynomial extrapolation.

Outputs:
- `c`: Linear model / filter
- `sums`: The extrapolated value applied to each window
- `resid`: The least squares residual
- `xs`: A time series `x[:,n] = Fⁿ(x[:,0])`
- `hs`: The observations `h[:,n] = h(x[:,n])`
- `rnorm`: The norm of resid
- `K`: The final degree of the filter
- `history`: Returned if `iterative=true`. The history of the final LSQR
  iteration
"""
function adaptive_birkhoff_extrapolation(h::Function, F::Function,
                    x0::AbstractVector; rtol::Number=1e-10, Kinit = 20,
                    Kmax = 100, Kstride=20, iterative::Bool=false,
                    Nfactor::Number=5, rre::Bool=true, ϵ::Number=0.0, weighted::Bool=false,
                    ortho::Bool=true)
    #
    d = length(h(x0));
    K = Kinit-Kstride
    N = ceil(Int, 2*Nfactor*K / d);

    sol, xs = nothing, nothing

    rnorm = Inf; 
    while (K+Kstride <= Kmax) && (rnorm > rtol)
        K += Kstride;
        N = ceil(Int, 2*Nfactor*K / d);

        sol = birkhoff_extrapolation(h, F, x0, N, K; iterative, x_prev=xs, rre, ϵ, weighted, ortho)
        rnorm = sol.resid_RRE
    end

    return sol
end

"""
    sum_stats(sums)

Given a sequence of sums applied to a filter (an output of invariant circle
extrapolation), find the average and standard deviation of the sums. Can be used
as a measure of how "good" the filter is.
"""
function sum_stats(sums)
    d, N = size(sums);
    sum_ave = sums*ones(N)/N;
    sum_res = zeros(d, N);
    for ii = 1:N
        sum_res[:, ii] = sums[:, ii] - sum_ave;
    end

    sum_cov = sum_res * sum_res' ./ (N-1)
    sum_std = sqrt.(eigvals(sum_cov));

    return sum_ave, sum_std
end

# This is a convolution. It could probably be done with FFTs
"""
    get_sum_ave(hs, c)

Given a signal `hs` (of size either d×N or N×1) and a filter `c` of size M×1,
compute the average of the M-N windows of the filter applied to the signal.
"""
function get_sum_ave(hs::AbstractArray, c::AbstractVector)
    d, N = size(hs);
    if N==1
        N = d;
        d = 1;
    end

    M = length(c);
    sum = d == 1 ? 0. : zeros(d);
    for ii = 1:N-M+1
        sum = sum + hs[:, ii:ii+M-1]*c;
    end
    sum = sum ./ (N-M+1);
    return sum
end

function find_closest(ω, ωs, not_assigned, tol)
    closest = 0;
    d_closest = 10;
    for ii = 1:length(ωs)
        if not_assigned[ii]
            d = abs(mod(ωs[ii] - ω + π, 2π) - π);
            replace = (d < d_closest) && (d < tol)
            d = replace ? d : d_closest;
            closest = replace ? ii : closest;
        end
    end

    return closest
end

function assign_modes(ω0, ωs, tol)
    N = length(ωs)
    p = zeros(Int, N);
    not_assigned = trues(N);

    for ii = 1:N÷2
        k1 = N÷2 + ii
        k2 = N - k1 + 1;

        jj = find_closest(ii*ω0, ωs, not_assigned, tol)
        if jj != 0
            not_assigned[jj] = false;
            p[k1] = jj;
        end

        jj = find_closest(-ii*ω0, ωs, not_assigned, tol)
        if jj != 0
            not_assigned[jj] = false;
            p[k2] = jj;
        end
    end

    return p[N÷2+1:end]
end

# Copied from pseudocode in Sander & Meiss
function SmallDenom(x, δ)
    x = mod(x, 1);
    n, d = 0, 1;
    pl, ql = 0, 1;
    pr, qr = 1, 0;

    while d==2 ? (abs(x-n/d) ≥ sqrt(δ)) : (abs(x - n/d) ≥ δ)
        n,d = pl+pr, ql+qr # Find the mediant
        if x < n/d
            pr, qr = n, d
        else
            pl, ql = n, d
        end
    end
    return (n, d);
end

function find_rationals(ωs::AbstractVector, dmax::Integer, tol::Number)
    Nisland = 1;
    ind = falses(length(ωs))
    for (ii, ω) in enumerate(ωs)
        n, d = SmallDenom(ω, tol)
        if d ≤ dmax
            Nisland = max(d, Nisland);
            ind[ii] = true;
        end
    end
    return ind, Nisland;
end

"""
    chebyshev_companion_matrix(v::AbstractVector)

Input:
- `v`: The coefficients of a Chebyshev polynomial

Output:
- `C`: The Chebyshev companion (colleague) matrix of polynomial `v`
"""
function chebyshev_companion_matrix(v::AbstractVector)
    K = length(v)-1
    C = zeros(K, K)

    if K == 0
        return C
    elseif K == 1
        C[1,1] = -v[1]/v[end]
        return C
    end

    C[1,2] = 1.
    C[2,1] = 0.5;
    for ii = 2:K-1
        C[ii, ii+1] = 0.5
        C[ii+1, ii] = 0.5
    end

    C[end, 1:end] += -v[1:end-1] ./ (2 .* v[end])
    C
end

"""
    palindromic_to_chebyshev(c::AbstractVector)

Input:
- `c`: A palindromic set of monomial coefficients for polynomial `p(z)`

Output:
- `v`: The coefficients of a Chebyshev polynomial `q` s.t.
    `q((z+inv(z))/2) = p(z)`
"""
function palindromic_to_chebyshev(c::AbstractVector)
    K = length(c)÷2
    @assert length(c) == 2K+1
    v = c[K+1:end]
    v[2:end] .*= 2.
    v
end

"""
    palindromic_cheb_roots(c::AbstractVector)

Get the roots of the chebyshev polynomial associated with the palindromic
polynomial with coefficients `c`.

Input:
- `c`: A palindromic set of monomial coefficients for polynomial `p(z)`

Output:
- `v`: The roots of the polynomial `q` satisfying `q((z+inv(z))/2) = p(z)`
"""
function palindromic_cheb_roots(c::AbstractVector{T}) where T
    v = palindromic_to_chebyshev(c)
    N = length(v)
    nnz = findlast((x) -> (x != zero(T)), v)
    v = v[1:nnz]

    C = chebyshev_companion_matrix(v)
    vcat(eigvals(C), Inf.*ones(N - nnz))
end


"""
    cheb_roots_to_roots(μs)

Input:
- `μs`: The roots of the polynomial satisfying `q((z+inv(z))/2) = p(z)`, where
   `p` is a palindromic polynomial

Output:
- `λs`: The roots of `p`
"""
function cheb_roots_to_roots(μs::AbstractVector)
    K = length(μs)
    λs = zeros(Complex{Float64}, 2K)
    for ii = 1:K
        μi = μs[ii]
        if μi == Inf
            λs[2ii-1:2ii] = [0.0, Inf]
        else
            λs[2ii-1:2ii] = roots(Polynomial([0.5, -μs[ii], 0.5]))
        end
    end

    λs
end

function eigenvalues_and_project(hs, c; growth_tol=1e-5)
    d, N = size(hs);
    sum_ave =  get_sum_ave(hs, c)
    λs = cheb_roots_to_roots(palindromic_cheb_roots(c));

    if (length(c) == 1) | (λs[1] == zero(typeof(λs[1]))) # Edge case: constant polynomial
        return [0.], [0.], sum_ave, zeros(1,2)
    end

    λs = λs[abs.(abs.(λs) .- 1) .< growth_tol]
    if length(λs) == 0
        @warn "No eigenvalues were found on the unit circle. This means rattol was probably too large"
    end
    K = length(λs);

    # eigenvectors
    vs = [λ^n for n in 0:N-1, λ in λs]

    # Getting the error (h minus the average)
    es = zeros(d, N);
    for ii = 1:d
        es[ii, :] = hs[ii,:] .- sum_ave[ii]
    end
    # Coefficient of the error
    coefs = vs \ es'

    # Norms of the coefficients
    norms = [norm(row) for row in eachrow(coefs)]
    p = sortperm(norms);
    p = p[end:-1:1]
    λs = λs[p]
    coefs = coefs[p, :]
    norms = norms[p]
    ωs = angle.(λs) ./ (2π);

    return λs, ωs, sum_ave, coefs
end


function get_w_grid(w0::AbstractVector, Ngrids::AbstractVector)
    d = length(w0) 

    if d < 1
        return [0.]
    elseif d == 1
        Ngrid = Ngrids[1]
        grid = [mid_mod(n*w0[1]) for n in -Ngrid:Ngrid]
        ind = sortperm(grid)

        return grid[ind]
    else
        grid_ld = get_w_grid(w0[1:d-1], Ngrids[1:d-1])
        grid_d  = get_w_grid(w0[d:d], Ngrids[d:d])

        grid = [mid_mod(x+y) for x in grid_ld, y in grid_d][:]
        ind = sortperm(grid)
        return grid[ind]
    end
end


function get_wk2_grid(w0::AbstractVector, Ngrids::AbstractVector)
    d = length(w0) 

    if d < 1
        return [0.], [1]
    elseif d == 1
        w = w0[1]
        Ngrid = Ngrids[1]
        if typeof(w) <: Rational
            Nisland = denominator(w)
            grid = [mid_mod(n*w) for n = 0:Nisland-1]
            ind = sortperm(grid)
            k2grid = zeros(Integer, 1, Nisland)
            return grid[ind], k2grid[ind]
        else
            grid = [mid_mod(n*w) for n in -Ngrid:Ngrid]; 
            ind = sortperm(grid)
            k2grid = zeros(Integer, 1,2Ngrid+1)
            k2grid[:] = (-Ngrid:Ngrid)[ind].^2
            return grid[ind], k2grid
        end
    elseif w0[1]==1//1
        return get_wk2_grid(w0[2:end], Ngrids[2:end])
    else
        grid_ld, k2grid_ld = get_wk2_grid(w0[1:d-1], Ngrids[1:d-1])
        grid_d , k2grid_d  = get_wk2_grid(w0[d:d], Ngrids[d:d])

        grid = [mid_mod(x + y) for x in grid_ld, y in grid_d][:]
        k2grid = vec([k1 + k2 for k1 in k2grid_ld, k2 in k2grid_d])

        ind = sortperm(grid)
        return grid[ind], k2grid[ind]
    end
end


function get_wk_grid(w0::AbstractVector, Ngrids::AbstractVector)
    d = length(w0) 

    if d < 1
        return [0.], [0]
    elseif d == 1
        w = w0[1]
        Ngrid = Ngrids[1]
        if typeof(w) <: Rational
            Nisland = denominator(w)
            grid = [mid_mod(n*w) for n = 0:Nisland-1]
            ind = sortperm(grid)
            kgrid = zeros(Integer, 1, Nisland)
            # kgrid[:] = (0:Nisland-1)[ind]
            return grid[ind], kgrid
        else
            grid = [mid_mod(n*w) for n in -Ngrid:Ngrid]; 
            ind = sortperm(grid)
            kgrid = zeros(Integer, 1,2Ngrid+1)
            kgrid[:] = (-Ngrid:Ngrid)[ind]
            return grid[ind], kgrid
        end
    else
        grid_ld, kgrid_ld = get_wk_grid(w0[1:d-1], Ngrids[1:d-1])
        grid_d , kgrid_d  = get_wk_grid(w0[d:d], Ngrids[d:d])

        grid = [mid_mod(x + y) for x in grid_ld, y in grid_d][:]

        Nk_ld = size(kgrid_ld,2)
        Nk_d = size(kgrid_d,2)
        kgrid = zeros(Integer, d, Nk_ld, Nk_d)

        @views for ii = 1:Nk_ld, jj = 1:Nk_d
            kgrid[1:d-1,ii,jj] = kgrid_ld[:,ii]
            kgrid[d,ii,jj] = kgrid_d[jj]
        end
        kgrid = reshape(kgrid,d,Nk_ld*Nk_d)
        
        ind = sortperm(grid)
        return grid[ind], kgrid[:,ind]
    end
end

function search_sorted_freqs(wgrid::AbstractVector, w::Number)
    N = length(wgrid)
    ind = searchsortedfirst(wgrid,w)

    k = 0
    if (ind == 1) || (ind == N+1)
        k = (abs(-abs(w)-wgrid[1]) > abs(abs(w)-wgrid[N])) ? N : 1
    elseif ind ≤ N
        k = (abs(w-wgrid[ind]) > abs(w-wgrid[ind-1])) ? ind-1 : ind
    end
    
    k
end

function match_ws(w0::AbstractVector, ws::AbstractVector, tol::Number, Ngrids::AbstractVector)

    wgrid = get_w_grid(w0, Ngrids)
    nearest = [search_sorted_freqs(wgrid, w) for w in ws]
    matches = [abs(mid_mod(ws[ii] - wgrid[nearest[ii]]))<tol for ii = 1:length(ws)]
    Nmatches = sum(matches)
    
    return matches, Nmatches
end

# Recursively finds a value for w0
# Returns the first frequency which exactly matches the vector `ws`
function initial_w0(w0::AbstractVector, ws::AbstractVector, Nw0::Integer, tol::Number, Ngrids::AbstractVector)
    w0_best = []
    matches_best = []
    Nmatches_best = 0  # Number of matches
    matches, Nmatches = match_ws(w0, ws, tol, Ngrids)

    if length(w0) == Nw0
        # If w0 is fully specified, return it and how well it did
        return w0, matches, Nmatches
    elseif length(w0) == 0
        # If none of w0 is specified yet, check if there is a rational frequency
    end

    # If w0 is partially specified, guess the next frequency in the vector
    for w in ws[.!matches]
        w0_candidate = vcat(w0, w)
        w0_candidate, matches, Nmatches = initial_w0(w0_candidate, ws, Nw0, tol, Ngrids)
        
        # If this candidate frequency specifies every other, use it
        if Nmatches == length(ws)
            return w0_candidate, matches, Nmatches
        end

        # Otherwise, if it matches more frequencies than the previous best, use it
        if Nmatches > Nmatches_best
            w0_best = w0_candidate
            matches_best = matches
            Nmatches_best = Nmatches
        end
    end

    return w0_best, matches_best, Nmatches_best
end

function sphere_volume(Nw0::Integer)
    if Nw0 == 1
        return 2
    elseif Nw0 == 2
        return π
    elseif Nw0 == 3
        return 4π/3
    elseif Nw0 == 4
        return (π^2)/2
    elseif Nw0 == 5
        return (π^2) * (8/15)
    elseif Nw0 == 6
        return (π^3)/6
    end

    @warn "Rotation Vector too long: need to implement better sphere volume routine"
    return -Inf
end

# Find approximate values for the smoothness by taking a best-fit line of the 
# Fourier coefficients.
function smoothness_estimate(h2norms::AbstractVector, Nw0::Integer, Nh::Integer, Nisland::Integer; keep_prop::Number=1/4)
    N = floor(Integer, length(h2norms)*keep_prop) 
    while (N >= 2) && (h2norms[2N]/h2norms[1])< 1e-10
        N = N-2
    end
    
    js = 0:N-1
    A = ones(N, 2)
    V = sphere_volume(Nw0)
    A[:,2] = -2 .* (js ./ (V*Nisland)) .^ (1/Nw0);
    
    v = (A) \ (log.(h2norms[1:N]))

    C_h = sqrt(exp(v[1]) / Nh)
    r_h = v[2]

    C_h, r_h
end

function wind_increment(wind, Nw, k)
    if (k==1) || (wind[k] < Nw)
        wind[k] = wind[k] + 1
        return wind[k]
    else
        wind[k] = wind_increment(wind, Nw-1, k-1) + 1
    end
end

# function chi2_pdf(x::Number, sigma::Number, k::Number)
#     den = 2^(k/2) * gamma(k/2);
#     x^(k/2-1)*exp(-x/2)/den
# end

# function chi2_cdf(x::Number, sigma::Number, k::Number)

# end

function w0_logposterior(w0::Number, w::Number, h2norm::Number, magk::Number, 
                         sigma_w::Number, r_h::Number, C_h::Number, Nh::Integer, h2min::Number)
    
    # As always, the w0=1/2 case is sensitive
    if w0 == 1/2
        sigma_w = sqrt(sigma_w)
    end

    # Probability of observing the frequency
    wdist = Normal()
    Pw = logpdf(wdist, mid_mod(w0-w)/sigma_w) - log(sigma_w)
    
    # Probability of observing the coefficient
    hdist = Chisq(Nh)
    sigma_h2 = (C_h*exp(-r_h*magk))^2
    Ph = logpdf(hdist, h2norm/sigma_h2) - log(sigma_h2)

    # Discount for not using other frequencies
    P_h2_lessthan_h2min = logcdf(hdist, h2min/sigma_h2)

    return Pw + Ph - P_h2_lessthan_h2min
end

function w0_logposterior(w0::AbstractVector, ws::AbstractVector, h2norms::AbstractVector, 
                         sigma_w::Number, C_h::Number, r_h::Number, N_h::Number, 
                         Ngrids::AbstractVector, searchwidth::Number)
    wgrid, k2grid = get_wk2_grid(w0, Ngrids)
    logposterior = 0.
    popped = falses(length(wgrid)) # Maybe this is inefficient, but whatev. It is 64x more efficient than wgrid

    Nw = length(wgrid)
    h2min = minimum(h2norms)
    
    for (ii,w) in enumerate(ws)
        jj_closest = search_sorted_freqs(wgrid, w) 
        jjs = mod1.((-searchwidth+jj_closest):(searchwidth+jj_closest), Nw)
        jjs = jjs[.!popped[jjs]]
        
        ## Needs to use the function with the correct arguments
        logposteriors = [w0_logposterior(wgrid[jj], w, h2norms[ii], sqrt(k2grid[jj]), sigma_w, 
                                         r_h, C_h, N_h, h2min) for jj in jjs]

        
        new_searchwidth=searchwidth
        while length(logposteriors) == 0 
            new_searchwidth = 2new_searchwidth
            jj_closest = search_sorted_freqs(wgrid, w) 
            jjs = mod1.((-new_searchwidth+jj_closest):(new_searchwidth+jj_closest), Nw)
            jjs = jjs[.!popped[jjs]]

            logposteriors = [w0_logposterior(wgrid[jj], w, h2norms[ii], sqrt(k2grid[jj]), sigma_w, 
                                         r_h, C_h, N_h, h2min) for jj in jjs]
        end

        kmax = argmax(logposteriors)
        logposterior = logposterior + logposteriors[kmax]
        popped[jjs[kmax]] = true
    end

    logposterior
end

function initial_w0(ws::AbstractVector, h2norms::AbstractVector, Nw0::Integer, Nisland::Integer, 
                    sigma_w::Number, r_h::Number, C_h::Number, N_h::Integer, Ngrids::AbstractVector; 
                    searchwidth = 5)
    Nw = length(ws)
    wind = [ii for ii in 1:Nw0]
    w0_best = vcat(Any[1//Nisland], zeros(Nw0))
    logposterior_best = -Inf
    finished = false
    
    # Loop through all possible values of w0 from the considered frequencies
    while !finished
        # Get the log posterior
        w0_candidate = vcat(Any[1//Nisland], ws[wind])
        logposterior = w0_logposterior(w0_candidate, ws, h2norms, sigma_w, r_h, C_h, N_h, 
                                       Ngrids,searchwidth)
        # If this is more likely, use it
        if logposterior > logposterior_best
            w0_best = w0_candidate
            logposterior_best=logposterior
        end

        if wind == vec(Nw-Nw0+1:Nw)
            finished=true
        end
        wind_increment(wind, Nw, Nw0)
    end

    return w0_best, logposterior_best
end

function refine_w0(w0::AbstractVector, ws::AbstractVector, coefs::AbstractArray, tol::Number, 
                   gridratio::Number)
    # Get wavenumber information
    Ngrids = floor(Integer,gridratio*sqrt(length(ws)/2)) * ones(Integer,length(w0))
    wgrid, kgrid = get_wk_grid(w0, Ngrids)
    nearest = [search_sorted_freqs(wgrid, w) for w in ws]
    matches = [abs(mid_mod(ws[ii]-wgrid[nearest[ii]]))<tol for ii = 1:length(ws)]
    
    coefs = coefs[matches]
    ks = kgrid[2:end, nearest[matches]] # Do not use island modenumbers in norm

    # Get loop minimizing basis
    hk = ks * Diagonal(coefs)
    H = real.(hk*hk')  # This works for islands too b/c of Parseval's Thm on the discrete frequency
    U = Matrix(cholesky(H).U)
    _, A = LLLplus.hkz(U)
    
    new_w0 = copy(w0)
    v = zeros(length(w0)-1)
    v[:] = w0[2:end]
    new_w0[2:end] = mid_mod.(A\v)
    D = Diagonal(Integer.(sign.(new_w0)));
    
    D*new_w0 # , D*A'*ks
end


"""
    get_w0(hs::AbstractArray, c::AbstractVector, Nw0::Number; matchtol::Number=1, 
                Nsearch::Integer=20, gridratio::Number=0., Nkz::Integer=50, sigma_w::Number = 1e-10)

Find the rotation vector from an output of Birkhoff RRE.

Input: 
- `hs`: The observable output from [`adaptive_birkhoff_extrapolation`](@ref)
- `c`: The filter output from [`adaptive_birkhoff_extrapolation`](@ref)
- `Nw0`: The dimension of the invariant torus
- `matchtol=1`: The tolerance at which we determine a frequency is an integer multiple of another 
   frequency for the refinement step. Can be increased/decreased depending on how well resolved `c` 
   is. For instance, if the torus is very complicated or the map is noisy, a larger tolerance may be 
   needed (say 1e-5).
- `(Nsearch,gridratio)=(20,0.)`: `Nsearch` is the number of independent roots of `c` (i.e. ±ω are the 
   same root) that we consider for choosing ω0. `gridratio` is gives the size of grid we use for the
   discrete optimization to find ω0. If the torus is dominated by higher harmonics, then one or both
   of these parameters may need to be increased. The default values (chosen by `gridratio==0.`)
   are `2.` for the 1D case and `1.` for greater than 1D.
- `Nkz=100`: Number of frequencies to use for finding the KZ basis. This value likely does not need
  to be tuned by the user.
- `sigma_w`: Frequency accuracy used in the Bayesian inference

Output:
- `w0`: The rotation vector. In the case of an invariant torus, is is of length `Nw0`. In the case 
  of an island, it is of length `Nw0+1`, where the first entry is rational and the rotation vector 
  is in elements `w0[2:Nw0+1]`.
- `logposterior`: The log posterior of the Bayesian determination of w0
"""
function get_w0(hs::AbstractArray, c::AbstractVector, Nw0::Number; matchtol::Number=1, 
                Nsearch::Integer=20, Nsearchcutoff::Number=1e-6, gridratio::Number=0., 
                Nkz::Integer=50, sigma_w::Number = 1e-10, rattol::Number = 1e-6, 
                maxNisland::Number = 10)
    # Get correct default grid ratio
    if gridratio==0.
        gridratio = (Nw0 == 1) ? 2. : 10.
    end

    ## Step 1: Get the roots sorted by the value of coefs
    _, ws, _, coefs = eigenvalues_and_project(hs, c)

    ## Step 1.5: Check if there is a rational rotation number
    H2norms = [real(H'*H) for H in eachrow(coefs)]
    while (Nsearch >= 2) && ((H2norms[2Nsearch]/H2norms[1]) < Nsearchcutoff) # Make sure that we aren't using garbage
        Nsearch = Nsearch-1
    end
    _, Nisland = find_rationals(ws[1:2:2Nsearch], maxNisland, rattol)
    
    ## Step 2: Get a candidate value of w0
    N_h = length(coefs[1,:])
    C_h, r_h = SymplecticMapTools.smoothness_estimate(H2norms, Nw0, N_h, Nisland)
    Ngrids = vcat([Nisland], floor(Integer, sqrt(2Nsearch)*gridratio) .* ones(Integer, Nw0))
    w0, logposterior = initial_w0(ws[1:2:2Nsearch], H2norms[1:2:2Nsearch], Nw0, Nisland, sigma_w, 
                                  r_h, C_h, N_h, Ngrids)

    ## Step 3: Find the optimal rotation vector by KZ basis
    Nkz = min(Nkz,size(coefs,1)÷4)
    coef_norms = [norm(row) for row in eachrow(coefs[1:2Nkz,:])]
    w0 = refine_w0(w0, ws[1:2Nkz], coef_norms, matchtol, gridratio)

    return w0, logposterior
end

function get_w0!(sol::BRREsolution, Nw0::Integer; matchtol::Number=1, 
    Nsearch::Integer=20, Nsearchcutoff::Number=1e-6, gridratio::Number=0., 
    Nkz::Integer=50, sigma_w::Number = 1e-10, rattol::Number = 1e-6, 
    maxNisland::Number = 10)
    w0, resid = get_w0(sol.hs, sol.c, Nw0; matchtol, Nsearch, gridratio, Nkz, sigma_w,
                       Nsearchcutoff,rattol,maxNisland)
    
    sol.d = Nw0
    sol.w0 = w0
    sol.resid_w0 = resid
end


mutable struct AdaptiveQR{T}
    N::Integer;
    M::Integer;
    D::Integer;
    N_val::Integer;
    factors::AbstractMatrix{T};
    τ::AbstractVector{T};
    X::AbstractMatrix{T};
    kinds::AbstractVector;
    A_val::AbstractMatrix;

    function AdaptiveQR(T::Type, N::Number, D::Number, N_val::Number)
        new{T}(N,
               0,
               D,
               N_val,
               zeros(T,N,1),
               zeros(T,1),
               zeros(T,1,D),
               zeros(Integer,1),
               zeros(T,N_val,1))
    end
end

function get_Q(QR::AdaptiveQR)
    @views LinearAlgebra.QR(QR.factors[:,1:QR.M], QR.τ[1:QR.M]).Q
end

function get_R(QR::AdaptiveQR)
    @views UpperTriangular(QR.factors[1:QR.M,1:QR.M])
end

function get_X(QR::AdaptiveQR)
    @views QR.X[1:QR.M,:]
end

function get_A_val(QR::AdaptiveQR)
    @views QR.A_val[:,1:QR.M]
end

function get_kinds(QR::AdaptiveQR)
    @views QR.kinds[1:QR.M]
end

function double!(QR::AdaptiveQR{T}) where T
    factors = QR.factors
    τ = QR.τ
    X = get_X(QR)
    kinds = get_kinds(QR)
    A_val = get_A_val(QR)

    M = QR.M
    N = QR.N
    D = QR.D
    N_val = QR.N_val

    buf = length(τ)


    QR.factors = Matrix{T}(undef, N, 2buf)
    QR.τ       = zeros(T, 2buf)
    QR.X       = Matrix{T}(undef, 2buf, D)
    QR.kinds   = Vector{Integer}(undef, 2buf)
    QR.A_val   = Matrix{T}(undef, N_val, 2buf)
    
    QR.factors[:,1:M] .= factors[:,1:M]
    QR.τ[1:M]         .= τ[1:M]
    QR.X[1:M, :]      .= X
    QR.kinds[1:M]     .= kinds
    QR.A_val[:,1:M]   .= A_val
end

function update_inds(ks_bank::AbstractArray, kinds::AbstractVector, sz::AbstractVector, dir::Number)
    k_iter = trues(size(ks_bank,2))
    k_iter[kinds] .= false
    sz_dir = copy(sz)
    sz_dir[dir] = sz_dir[dir] + 1

    for (ii, k) in enumerate(eachcol(ks_bank))
        if k_iter[ii]
            if !(prod(abs.(k) .<= sz_dir))
                k_iter[ii] = false
            end
        end
    end

    kinds_out = (1:size(ks_bank,2))[k_iter]
    kinds_out
end

function update!(QR::AdaptiveQR, A::AbstractArray, kinds::AbstractArray, A_val::AbstractArray)
    M = QR.M
    N = QR.N
    M2 = size(A,2)
    @assert size(A,1) == N

    while QR.M+M2 > length(QR.τ)
        double!(QR)
    end

    ind1 = 1:M
    not_ind1 = M+1:N;
    ind2 = M+1:M+M2
    
    Q = get_Q(QR)
    Q_proj = Q'*A

    QR.factors[ind1, ind2] .= Q_proj[ind1,:];
    factor_up, τ_up = LAPACK.geqrf!(Q_proj[not_ind1,:])
    
    QR.factors[not_ind1, ind2] .= factor_up
    QR.τ[ind2] .= τ_up

    QR.kinds[ind2] .= kinds

    QR.A_val[:,ind2] .= A_val

    QR.M = M+M2;

end

function downdate!(QR::AdaptiveQR, M2::Number)
    QR.M = QR.M-M2;
end

function solve!(QR::AdaptiveQR{T}, B::AbstractArray{T}) where T
    R = get_R(QR)
    Q = get_Q(QR)
    M = QR.M

    @views QR.X[1:M, :] .= R\(Q'*B)[1:M,:]
end


function get_update_matrix(w0::AbstractVector, N::Number, n::Number, sz::AbstractVector)
    Nw0 = length(w0)
    if Nw0 == 1
        M = (n==1) ? 2 : 2sz[1]+1
        A = zeros(ComplexF64, N, M)
        ks = zeros(Integer, 1, M)

        if n == 1
            # lambdainv = exp.(-2π*im*w0[1]*(sz[1]+1))
            # lambda = exp.(2π*im*w0[1]*(sz[1]+1))
            f = 2π*im*w0[1]*(sz[1]+1)

            A[:,1] .= exp.(-f .* (0:N-1))# [lambdainv^t for t = 0:N-1]
            # A[:,2] = [lambda^t for t = 0:N-1]
            A[:,2] .= conj.(A[:,1])

            ks[:] .= [-(sz[1]+1),sz[1]+1]
        else
            for (ii,k) in enumerate(-sz[1]:sz[1])
                f = 2π*im*w0[1]*k
                A[:,ii] .= exp.(f .* (0:N-1))
                ks[ii] = k;
            end
        end
        
        return A, ks
    end

    A1, ks1 = get_update_matrix(w0[1:end-1], N, n, sz[1:end-1])
    A2, ks2 = get_update_matrix(w0[[Nw0]], N, n-Nw0+1, sz[[Nw0]])

    M1 = size(A1,2)
    M2 = size(A2,2)
    M = M1*M2

    # A = zeros(ComplexF64, N, M)
    A = Matrix{ComplexF64}(undef, N, M)
    ks = zeros(Integer, Nw0, M)

    for ii = 1:M2
        ind = (ii-1)*M1 .+ (1:M1)
        A[:,ind] .= Diagonal(A2[:,ii]) * A1 
        ks[1:end-1,ind] .= ks1
        ks[end,ind] .= ks2[ii]
    end
    
    A, ks
end

function get_torus_err(B_val::AbstractArray, QR::AdaptiveQR, eps_ada::Number, ks_bank::AbstractArray)
    X = get_X(QR)
    A_val = get_A_val(QR)
    
    v = A_val * X - B_val
    knorms = [sqrt(col'col + 1) for col in eachcol(@view ks_bank[:,get_kinds(QR)])]
    sqrt(norm(v)^2 + eps_ada*norm(Diagonal(knorms)*X)^2)
end


# """
    
# """
function adaptive_get_torus(w0::AbstractVector, hs::AbstractMatrix; 
                            validation_fraction::Number=0.05, ridge_factor::Number=10,
                            eps_ada::Number=1e-8)
    
    # Preprocess for an island
    Nisland = 1
    Nw0 = length(w0)
    if typeof(w0[1]) <: Rational
        Nisland = denominator(w0[1])
        w0_old = w0
        w0 = zeros(Nw0-1)
        w0[:] = abs.(mid_mod.(w0_old[2:end] .* Nisland))
        Nw0 = Nw0-1
    end

    D_pre, N_pre = size(hs)
    hs = reshape(hs[:,1:(N_pre÷Nisland) * Nisland], D_pre*Nisland, N_pre÷Nisland) # Reshape to island
    D, N = size(hs)
    
    # Problem dimensions
    
    N_train = floor(Integer, (1-validation_fraction)*N)

    ind_train = 1:N_train
    ind_val = N_train+1:N

    # Initialize the least-squares problem
    sz = zeros(Integer, Nw0)
    B = hs'
    B_train = B[ind_train, :] .* (1. + 0. * im)
    B_val   = B[ind_val, :]   .* (1. + 0. * im)
    err_den = norm(B[ind_val, :]) # normalization for the relative error
    
    mode_bank = ones(ComplexF64, N, 1)
    ks_bank = zeros(Integer, Nw0, 1)
    for ii = 1:Nw0
        sz_ii = zeros(Integer, Nw0)
        sz_ii[1:ii-1] .= 1;
        Aii, ksii = get_update_matrix(w0,N,ii,sz_ii)

        mode_bank = hcat(mode_bank, Aii)
        ks_bank   = hcat(ks_bank,  ksii)
    end
    ind_A=[1]
    
    QRs = [AdaptiveQR(ComplexF64, N_train, D, N-N_train) for _ = 1:Nw0]
    for ii = 1:Nw0; 
        update!(QRs[ii], mode_bank[ind_train,ind_A],ind_A, mode_bank[ind_val, ind_A])
        solve!(QRs[ii],B_train)
    end

    dir_best = 1;
    err_best = get_torus_err(B_val, QRs[dir_best], eps_ada, ks_bank) / err_den
    M_best   = 1;

    Lmode = size(mode_bank,2);
    
    

    # Adaptively find the best validation error
    for _ = 1:N_train # This loop should always break, but making it finite for error catching
        # Consider updating in each dimension
        err_candidates = zeros(Nw0)

        for jj = 1:Nw0
            kind_update = update_inds((@view ks_bank[:,1:Lmode]), get_kinds(QRs[jj]), sz, jj)
            
            if length(kind_update)+QRs[jj].M > N_train
                err_candidates[jj] = Inf
            else
                @views update!(QRs[jj],(mode_bank[ind_train,kind_update]), kind_update, mode_bank[ind_val, kind_update])
                solve!(QRs[jj], B_train)
                err_candidates[jj] = get_torus_err(B_val, QRs[jj], eps_ada, ks_bank)/err_den
            end
        end
        
        # Update in the best direction or break
        dir = argmin(err_candidates)
        err_dir = err_candidates[dir]
        if err_dir < ridge_factor*err_best
            # If this is the best we've seen, update best values
            if err_dir < err_best
                err_best = err_dir
                dir_best = dir
                M_best   = QRs[dir].M
            end

            # Either way, we take a next stpe
            mode_next, ks_next = get_update_matrix(w0,N,dir,sz .+ 1)

            Lnext = size(mode_next,2)
            while Lnext+Lmode >= size(mode_bank,2)
                tmp = mode_bank
                tmp_k = ks_bank
                
                mode_bank = Matrix{ComplexF64}(undef, N, 2*size(tmp,2))
                ks_bank = zeros(Integer,Nw0,2*size(tmp,2))
                mode_bank[:,1:Lmode] .= tmp[:,1:Lmode]
                ks_bank[:,1:Lmode] .= tmp_k[:,1:Lmode]
            end
            mode_bank[:,Lmode+1:Lmode+Lnext] .= mode_next
            ks_bank[:,Lmode+1:Lmode+Lnext] .= ks_next

            Lmode = Lmode+Lnext
            sz[dir] = sz[dir]+1
        else
            break
        end
    end
    
    # Recover the best solution
    QR = QRs[dir_best]
    downdate!(QR,QR.M-M_best)
    solve!(QR,B_train)
    X = get_X(QR)
    kinds = get_kinds(QR)

    # Create the torus
    ks = ks_bank[:,kinds]
    Nks = [maximum(row) for row in eachrow(ks)]
    tor = FourierTorus(D_pre, Nks; τ = 2π .* w0, p=Nisland);
    for (ii, k) in enumerate(eachcol(ks))
        tor.a[:, :, (k + Nks .+ 1)...] = X[ii, :]
    end

    return tor, err_best
end

# max_size caps the size of linear system we are considering
function adaptive_get_torus!(sol::BRREsolution; 
                             max_size::Number=2000,
                             validation_fraction::Number=0.05, ridge_factor::Number=10,
                             eps_ada::Number=1e-8)
    h_sz_2 = Integer(min(max_size, size(sol.hs,2)))
    tor, resid_tor = adaptive_get_torus(sol.w0, sol.hs[:,1:h_sz_2]; 
                                        validation_fraction, ridge_factor, eps_ada)

    sol.tor = tor
    sol.resid_tor = resid_tor
end


"""
    get_torus(w0::AbstractVector, hs::AbstractArray, ws::AbstractVector, coefs::AbstractArray;
                   coef_cutoff::Number = 1e-6, tol::Number=1e-9, gridratio::Number=10., 
                   validation_fraction::Number=0.05)

Project a sequence of observables onto an `InvariantTorus`. To handle potential scale separation 
between dimensions (long thin tori), the outputs ws and coefs

Input:
- `w0`: The rotation vector, likely obtained via [`get_w0`](@ref)
- `hs`: The sequence of observables, likely taken from [`adaptive_birkhoff_extrapolation`](@ref)
- `ws`: The rotation number roots returned from [`get_w0`](@ref)
- `coefs`: The projected inner product coefficients from [`get_w0`](@ref)
- `coef_cutoff=1e-6`
"""
function get_torus(w0::AbstractVector, hs::AbstractArray, ws::AbstractVector, coefs::AbstractArray;
                   coef_cutoff::Number = 1e-6, tol::Number=1e-9, gridratio::Number=10., 
                   validation_fraction::Number=0.05)
    N = size(hs, 2)
    
    ## Step 1: Estimate the maximum number of modes in each direction
    coef_norms = [norm(x) for x in eachrow(coefs)]
    abs_cutoff = coef_cutoff * (norm(coef_norms) / sqrt(length(coef_norms)))
    N_est = sum(coef_norms .> abs_cutoff)
    N_est = 2*(N_est÷2)
    
    Ngrids = floor(Integer,gridratio*sqrt(N_est/2)) * ones(Integer,length(w0))
    wgrid, kgrid = get_wk_grid(w0, Ngrids)
    nearest = [search_sorted_freqs(wgrid, w) for w in ws[1:N_est]]
    matches = [abs(mid_mod(ws[ii]-wgrid[nearest[ii]]))<tol for ii = 1:N_est]
    
    ks = kgrid[:, nearest[matches]]
    Nks = [maximum(abs.(row)) for row in eachrow(ks)]

    ## Step 2: Project `hs` onto the modes
    
    # Get modenumbers
    ks = zeros(Integer, length(Nks),prod(2 .* Nks .+ 1))
    cart = CartesianIndices(Tuple([-Nk:Nk for Nk in Nks]))
    for ii = 1:length(Nks), jj in 1:prod(2 .* Nks .+ 1)
        ks[:,jj].=Tuple(cart[jj])
    end

    # Get modes
    fs = (2π*im) .* (ks'*w0)
    modes = [exp(f*t) for t in 0:N-1, f in fs]

    # Project onto training data and get validation error
    N_val = ceil(Integer, validation_fraction*N)
    
    ind_train = 1:N-N_val
    ind_val = N-N_val+1:N
    
    coefs = modes[ind_train,:]\hs[:, ind_train]'
    err_val = norm(modes[ind_val,:]*coefs - hs[:,ind_val]')/norm(modes[ind_val,:])


    ## Step 3: Form the invariant torus
    tor = FourierTorus(size(hs,1), Nks; τ = 2π .* w0);
    for (ii, k) in enumerate(eachcol(ks))
        tor.a[:, 1, (k + Nks .+ 1)...] = coefs[ii, :]
    end

    return tor, err_val
end

# OLD (probably)
function get_Nmode(ω0::Number, N::Number, modetol::Number)
    W = wba_weight(1,N);
    w = diag(W*W);

    V = [exp(im*m*n*ω0) for m in 1:N, n in 1:N];
    tmp = abs.(V*w)
    vals = cumsum(tmp);

    Nmode = findfirst((x)->(x>modetol), vals)
    Nmode = (2*(Nmode÷2) == Nmode) ? Nmode-1 : Nmode
    Nmodemax = (2*(N÷2) == N) ? N-1 : N

    return min(Nmode, Nmodemax)
end


function get_circle_coef(hs::AbstractArray, ω0::Number;
                         modetol::Number = 0.5, Nmode::Integer=-1)
    den = denoms(ContFrac(ω0/2π));
    d, N = size(hs);

    if Nmode == -1
        Nmode = get_Nmode(ω0, N, modetol)
    end
    
    # λ = exp(2π*im * ω0);
    ωs = zeros(Nmode);
    ωs[1] = 0;
    for ii = 1:Nmode÷2
        ωs[2ii] = ii*ω0
        ωs[2ii+1] = -ii*ω0
    end

    w = wba_weight(1, N)
    # w = ones(N)
    # w = w / sqrt(sum(w.^2))

    vs = Diagonal(w)*[exp(m*im*ωi) for m=0:N-1, ωi = ωs]
    rhs = Diagonal(w)*hs[:,1:N]'

    return vs \ rhs
end


function coefs_to_circle(sum_ave::AbstractArray, circle_coef::AbstractArray,
                         ω0::Number, Nisland::Integer)
    z = FourierCircle(size(circle_coef,1)÷2; τ=ω0, p=Nisland)

    for kk = 1:Nisland
        ind_kk = 2kk-1:2kk
        set_a0!(z, sum_ave[ind_kk]; i_circle=kk)

        for jj = 1:get_Na(z)
            # jj = p2[ii]
            Am = zeros(2,2);
            Am[:, 1] = 2 .* real.(circle_coef[2jj, ind_kk])
            Am[:, 2] = -2 .* imag.(circle_coef[2jj, ind_kk])
            set_Am!(z, jj, Am; i_circle=kk)
        end
    end

    z
end

"""
    get_circle_info(hs::AbstractArray, c::AbstractArray; rattol::Number=1e-8,
                    ratcutoff::Number=1e-4, max_island_d::Integer=30)

Get a Fourier representation of an invariant circle from the observations
`hs` and the learned filter `c` (see `adaptive_invariant_circle_model` and
`invariant_circle_model` to find the filter). See `get_circle_residual` for an
a posteriori validation of the circle.

Optional Arguments:
- `rattol`: Roots are judged to be rational if |ω-m/n|<rattol
- `ratcutoff`: Relative prominence needed by a linear mode to qualify as
  "important" for deciding whether the sequence is an island
- `max_island_d`: Maximum denominator considered for islands.
- `modetol`: Tolerance (typically less than 1) used for determining the number
  of Fourier modes. Higher numbers result in more Fourier modes at the expense
  of robustness of the least-squares system.  Decrease it (e.g. to 0.5) for
  noisy data

Output:
- `z`: An invariant circle of type `FourierCircle`
"""
function get_circle_info(hs::AbstractArray, c::AbstractArray;
                         rattol::Number=1e-8, ratcutoff::Number=1e-2,
                         max_island_d::Integer=30, ϵ::Number=0.0,
                         modetol::Number=0.5)
    λs, ωs, sum_ave, coefs = eigenvalues_and_project(hs, c)

    if max_island_d == 0
        K = length(c)÷2;
        max_island_d = K
    end
    # Check for a rational part

    coefnorm = norm(coefs);
    # Only check the modes for rational numbers that account for a significant
    # amount of the energy of the signal
    ind = [norm(coefs[ii, :])/coefnorm > ratcutoff for ii = 1:length(λs)]
    ind_rat, Nisland = find_rationals(ωs[ind], max_island_d, rattol)

    # Get the coefficients of the invariant circle
    ω0 = ωs[1]*2π
    circle_coef = [];
    if Nisland == 1
        # Find coefficients of circles
        circle_coef = get_circle_coef(hs, ω0; modetol)
    else
        # If there are islands, rerun extrapolation to find reduced filter
        d, N = size(hs);
        hs_island = reshape(hs[:, 1:(N÷Nisland)*Nisland], d*Nisland, N÷Nisland)
        zs = Vector{Any}(undef, Nisland);
        K = length(c) ÷ 2;
        Kisland = K ÷ Nisland +1;
        c_island, ~, resid = vector_rre_backslash(hs_island, Kisland; ϵ);

        λs, ωs, sum_ave, coefs = eigenvalues_and_project(hs_island, c_island)
        # Then, find coefficients of each island
        ω0 = ωs[1]*2π
        circle_coef = get_circle_coef(hs_island, ω0; modetol)
    end

    return coefs_to_circle(sum_ave, circle_coef, ω0, Nisland)
end
