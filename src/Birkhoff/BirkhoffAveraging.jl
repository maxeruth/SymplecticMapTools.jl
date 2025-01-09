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


"""
    birkhoff_extrapolation(h::Function, F::Function, x0::AbstractVector,
                           N::Integer, K::Integer; iterative::Bool=false,
                           x_prev::Union{AbstractArray,Nothing}=nothing,
                           rre::Bool=false)

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
                                rre::Bool=false, ϵ::Number=0.0)
    x = deepcopy(x0);
    h0 = h(x);
    d = length(h0);

    @assert N*d ≥ K

    Nx = N+2K+1
    hs = zeros(d, Nx);
    xs = zeros(d, Nx);

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
            c, sums, resid, history = vector_rre_iterative(hs, K; ϵ)
        else
            c, sums, resid = vector_rre_backslash(hs, K; ϵ)
        end
    else
        if iterative
            c, sums, resid, history = vector_mpe_iterative(hs, K; ϵ)
        else
            c, sums, resid = vector_mpe_backslash(hs, K; ϵ)
        end
    end

    return c, reshape(sums, d, N+1), reshape(resid, d, N), xs, hs, history;
end

function unorm(hs::AbstractArray)
    us = diff(hs, dims=2)
    us2 = [ui'*ui for ui in eachcol(us)]
    sqrt(weighted_birkhoff_average(us2))
end

"""
    adaptive_birkhoff_extrapolation(h::Function, F::Function,
                    x0::AbstractVector; rtol::Number=1e-10, Kinit = 20,
                    Kmax = 100, Kstride=20, iterative::Bool=false,
                    Nfactor::Integer=1, rre::Bool=true)

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
                    Nfactor::Number=1, rre::Bool=true, ϵ::Number=0.0)
    #
    d = length(h(x0));
    K = Kinit-Kstride
    N = ceil(Int, 2*Nfactor*K / d);

    c, sums, resid, xs, hs, history = 0,0,Inf,nothing,nothing,0;
    rnorm = Inf;
    while (K+Kstride <= Kmax) && (rnorm > rtol)
        K += Kstride;
        N = ceil(Int, 2*Nfactor*K / d);

        c, sums, resid, xs, hs, history = birkhoff_extrapolation(
                                      h, F, x0, N, K; iterative, x_prev=xs, rre, ϵ)
        rnorm = norm(resid)/unorm(hs)
    end

    return c, sums, resid, xs, hs, rnorm, K, history
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
        Ngrid = Ngrids[1]
        grid = [mid_mod(n*w0[1]) for n in -Ngrid:Ngrid]; # mod.((-Ngrid:Ngrid).*w0[1] .+ 0.5, 1) .- 0.5
        ind = sortperm(grid)
        k2grid = zeros(Integer, 2Ngrid+1)
        k2grid[:] = (-Ngrid:Ngrid)[ind]
        return grid[ind], kgrid
    else
        grid_ld, k2grid_ld = get_wk_grid(w0[1:d-1], Ngrids[1:d-1])
        grid_d , k2grid_d  = get_wk_grid(w0[d:d], Ngrids[d:d])

        grid = [mid_mod(x + y) for x in grid_ld, y in grid_d][:]
        k2grid = vec([k1^2 + k2^2 for k1 in k2grid_ld, k2 in k2grid_d])

        ind = sortperm(grid)
        return grid[ind], k2grid[ind]
    end
end


function get_wk_grid(w0::AbstractVector, Ngrids::AbstractVector)
    d = length(w0) 

    if d < 1
        return [0.], [1]
    elseif d == 1
        Ngrid = Ngrids[1]
        grid = [mid_mod(n*w0[1]) for n in -Ngrid:Ngrid]; #mod.((-Ngrid:Ngrid).*w0[1] .+ 0.5, 1) .- 0.5
        ind = sortperm(grid)
        kgrid = zeros(Integer, 1,2Ngrid+1)
        kgrid[:] = (-Ngrid:Ngrid)[ind]
        return grid[ind], kgrid
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
        # display(Ns)
        
        # ind = sum(abs.(kgrid),dims=1) .< maximum(Ngrids)
        # ind = reshape(ind,length(ind))
        # grid = grid[ind]
        # kgrid = kgrid[:,ind]

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

function var_h_normalization(r_h::Number, Ngrids::AbstractVector, Nw0::Number)
    w0 = randn(length(Ngrids)) # whatever
    _, k2grid = get_wk2_grid(w0, Ngrids::AbstractVector)
    magkgrid = sqrt.(k2grid)

    vars_h = exp.(-(2r_h) .* magkgrid)
    
    sum(vars_h)*Nw0
end

function wind_increment(wind, Nw, k)
    if (k==1) || (wind[k] < Nw)
        wind[k] = wind[k] + 1
        return wind[k]
    else
        wind[k] = wind_increment(wind, Nw-1, k-1) + 1
    end
end

function chi2_pdf(x::Number, sigma::Number, k::Number)
    den = 2^(k/2) * gamma(k/2);
    x^(k/2-1)*exp(-x/2)/den
end

function chi2_cdf(x::Number, sigma::Number, k::Number)

end

function w0_likelihood(w0::AbstractVector, w::Number, h2norm::Number, magk::Number, 
                       sigma_w::Number, r_h::Number, N_h::Number, h2min::Number)
    
    # Probability of observing the frequency
    wdist = Normal()
    Pw = logpdf(mid_mod(w0-w)/sigma_w)

    # Probability of observing the coefficient
    D = length(h)
    magh = norm(h)
    sigma_h = exp(-r_h*magk)/N_h
    Ph = -(d/2)*log(2π*sigma_h^2)

    # Probability of the coefficient being larger than or equal to the minimum coefficient
    hdist = Chisq(length(w0))
end

function w0_likelihood(w0::AbstractVector, ws::AbstractVector, h2norms::AbstractVector, sigma_w::Number, 
                       r_h::Number, N_h::Number, Ngrids::AbstractVector,searchwidth::Number)
    wgrid, k2grid = get_wk2_grid(w0, Ngrids)
    loglikelihood = 0.

    Nw = length(ws)
    h2min = minimum(h2norms)
    
    for w in ws
        jj_closest = search_sorted_freqs(wgrid, w) 
        jjs = mod1.(-width+jj_closest:width+jj_closest, Nw)
        
        loglikelihoods = [w0_likelihood(w0, wgrid[jj], h2norms[jj], sqrt(k2grid[jj]), sigma_w, r_h, 
                                        N_h, h2min) for jj in jjs]
        kmax = argmax(loglikelihoods)
        loglikelihood = loglikelihood + popat!(loglikelihoods, kmax)

    end
end

# Recursively finds a value for w0
# Returns the first frequency which exactly matches the vector `ws`
function initial_w0(ws::AbstractVector, hs::AbstractVector, Nw0::Integer, sigma_w::Number, 
                    r_h::Number, Ngrids::AbstractVector; searchwidth = 5)
    Nw = length(ws)
    wind = [ii for ii in 1:Nw0]
    w0_best = zeros(Nw0)
    loglikelihood_best = -Inf
    finished = false
    
    N_h = var_h_normalization(r_h, Ngrids, Nw0)
    h2norms = [h'*h for h in eachcol(hs)]
    # If w0 is partially specified, guess the next frequency in the vector
    while !finished
        w0_candidate = ws[wind]
        loglikelihood = w0_likelihood(w0_candidate, ws, h2norms, sigma_w, r_h, N_h, Ngrids,searchwidth)
        
        # If this is more likely, use it
        if loglikelihood > loglikelihood_best
            w0_best = w0_candidate
            loglikelihood_best=loglikelihood
        end

        if wind == vec(Nw-Nw0+1:Nw)
            finished=true
        end
        wind_increment(wind, Nw, Nw0)
    end

    return w0_best# , matches_best, Nmatches_best
end

function refine_w0(w0::AbstractVector, ws::AbstractVector, coefs::AbstractArray, tol::Number, 
                   gridratio::Number)
    # Get wavenumber information
    Ngrids = floor(Integer,gridratio*sqrt(length(ws)/2)) * ones(Integer,length(w0))
    wgrid, kgrid = get_wk_grid(w0, Ngrids)
    nearest = [search_sorted_freqs(wgrid, w) for w in ws]
    matches = [abs(mid_mod(ws[ii]-wgrid[nearest[ii]]))<tol for ii = 1:length(ws)]
    
    coefs = coefs[matches]
    ks = kgrid[:, nearest[matches]]

    # Get loop minimizing basis
    hk = ks * Diagonal(coefs)
    H = real.(hk*hk')
    U = Matrix(cholesky(H).U)
    _, A = LLLplus.hkz(U)
    
    new_w0 = A\w0
    D = Diagonal(Integer.(sign.(new_w0)));
    
    display(H)
    display(A'*H*A)

    D*new_w0 # , D*A'*ks
end


"""
    get_w0(hs::AbstractArray, c::AbstractVector, dim::Number; tol::Number=1e-9, 
                Nsearch::Integer=20, gridratio::Number=0., Nkz::Integer=50)

Find the rotation vector from an output of Birkhoff RRE.

Input: 
- `hs`: The observable output from [`adaptive_birkhoff_extrapolation`](@ref)
- `c`: The filter output from [`adaptive_birkhoff_extrapolation`](@ref)
- `dim`: The dimension of the invariant torus
- `tol=1e-8`: The tolerance at which we determine a frequency is an integer multiple of another 
   frequency. Can be increased/decreased depending on how well resolved `c` is. For instance, if 
   the torus is very complicated or the map is noisy, a larger tolerance may be needed (say 1e-5).
- `(Nsearch,gridratio)=(20,0.)`: `Nsearch` is the number of independent roots of `c` (i.e. ±ω are the 
   same root) that we consider for choosing ω0. `gridratio` is gives the size of grid we use for the
   discrete optimization to find ω0. If the torus is dominated by higher harmonics, then one or both
   of these parameters may need to be increased. The default values (chosen by `gridratio==0.`)
   are `2.` for the 1D case and `1.` for greater than 1D.
- `Nkz=100`: Number of frequencies to use for finding the KZ basis. This value likely does not need
  to be tuned by the user.

Output:
- `w0`: The rotation vector. In the case of an invariant torus, is is of length `dim`. In the case 
  of an island, it is of length `dim+1`, where the first entry is rational and the rotation vector 
  is in elements `w0[2:dim+1]`.
- `matchedratio`: The number of frequencies that were matched to tolerance `tol` over the total 
  number of considered frequencies (`Nsearch`). If it is `1.`, then all frequencies were matched.
  Otherwise, the number will be below `1.`. `matchedratio==1.` does not necessarily imply that the 
  rotation vector is correct. Errors could still include: `tol` is too large, `Nsearch` is too small 
  to see a basis of w0, `dim` is too large. However, this can be used as to check the reliability of 
  this individual step, while a final check of the torus parameterization is the final necessary 
  component of the algorithm.
- `ws`: The full set of frequencies associated with the roots of `c`. Both `ws` and `coefs` (below)
  are used as inputs to the function [`get_torus`](@ref).
- `coefs`: The inner product of the frequencies `ws` against the signal `hs`
"""
function get_w0(hs::AbstractArray, c::AbstractVector, dim::Number; tol::Number=1e-9, 
                Nsearch::Integer=20, gridratio::Number=0., Nkz::Integer=50)
    # Get correct default grid ratio
    if gridratio==0.
        gridratio = (dim == 1) ? 2. : 10.
    end

    # Step 0: Normalize the input
    NW = size(hs,2)
    W = wba_weight(1,NW)
    hnorm = sqrt(sum((hs.*hs*(W*W))[:]))
    hs = hs ./ hnorm;
    display("norm of hs = $(norm(hs)), NW = $(NW)")
    
    # Step 1: Get the roots sorted by the value of coefs
    _, ws, _, coefs = eigenvalues_and_project(hs, c)

    # Step 2: Get a candidate value of w0
    Ngrids = floor(Integer, sqrt(Nsearch)*gridratio) .* ones(Integer, dim)
    w0, matches, Nmatches = initial_w0(Number[], ws[1:2:2Nsearch], dim, tol, Ngrids)
    
    # Step 3: Find the optimal rotation vector by KZ basis
    coef_norms = [norm(row) for row in eachrow(coefs[1:2Nkz,:])]
    w0 = refine_w0(w0, ws[1:2Nkz], coef_norms, tol, gridratio)

    return w0, Nmatches/Nsearch, ws, coefs
end


mutable struct AdaptiveQR{T}
    N::Integer;
    M::Integer;
    D::Integer;
    Q::AbstractMatrix{T};
    R::AbstractMatrix{T};
    X::AbstractMatrix{T};
    kinds::AbstractVector;

    function AdaptiveQR(T::Type, N::Number, D::Number)
        new{T}(N,
               0,
               D,
               zeros(T,N,1),
               zeros(T,1,1),
               zeros(T,1,D),
               zeros(Integer,1))
    end
end

function get_Q(QR::AdaptiveQR)
    @views QR.Q[:,1:QR.M]
end

function get_R(QR::AdaptiveQR)
    @views UpperTriangular(QR.R[1:QR.M,1:QR.M])
end

function get_X(QR::AdaptiveQR)
    @views QR.X[1:QR.M,:]
end

function get_kinds(QR::AdaptiveQR)
    @views QR.kinds[1:QR.M]
end

function double!(QR::AdaptiveQR{T}) where T
    Q = get_Q(QR)
    R = get_R(QR)
    X = get_X(QR)
    kinds = get_kinds(QR)

    M = QR.M
    N = QR.N
    D = QR.D

    buf = size(QR.Q,2)

    QR.Q =  Matrix{T}(undef, N,    2buf)
    QR.R =  zeros( T,        2buf, 2buf)
    QR.X =  Matrix{T}(undef, 2buf, D)
    QR.kinds = Vector{Integer}(undef, 2buf)
    
    QR.Q[:,1:M]   .= Q
    QR.R[1:M,1:M] .= R
    QR.X[1:M, :]  .= X
    QR.kinds[1:M] .= kinds
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

function update!(QR::AdaptiveQR, A::AbstractArray, kinds::AbstractArray)
    M = QR.M
    N = QR.N
    M2 = size(A,2)
    @assert size(A,1) == N

    while QR.M+M2 > size(QR.R,2)
        double!(QR)
    end

    ind1 = 1:M
    ind2 = M+1:M+M2
    Q = get_Q(QR)
    QR.R[ind1, ind2] = Q'*A;
    A = A - Q * QR.R[ind1, ind2] 
    
    Q2, R2 = qr(A)
    QR.R[ind2,ind2] .= R2
    QR.Q[:,ind2] .= Matrix(Q2)
    QR.kinds[ind2] .= kinds

    QR.M = M+M2;
end

function downdate!(QR::AdaptiveQR, M2::Number)
    QR.M = QR.M-M2;
end

# function solve(QR::AdaptiveQR{T}, B::AbstractArray{T}) where T
#     # tmp = zeros(ComplexF64, QR.M,size(B,2))
#     # println("Bsz = $(size(get_Q(QR),2)*size(B,2))")
#     Q = get_Q(QR)
#     tmp = Q'*B
#     get_R(QR)\tmp
# end

# function solve!(X::AbstractArray,tmp::AbstractArray,QR::AdaptiveQR,B::AbstractArray)
#     # println("Bsz = $(size(get_Q(QR),2)*size(B,2))")
#     M = QR.M
    
#     tmp2 = @view tmp[1:M,:]
#     Q = get_Q(QR)
#     # tmp2[:,:] = Q'*B
#     mul!(tmp2,Q',B)
#     X[1:M,:] = get_R(QR)\tmp2
# end

function solve!(QR::AdaptiveQR{T}, B::AbstractArray{T}) where T
    Q = get_Q(QR)
    R = get_R(QR)
    M = QR.M
    
    QR.X[1:M, :] .= R\(Q'*B)
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

    A = zeros(ComplexF64, N, M)
    ks = zeros(Integer, Nw0, M)

    for ii = 1:M2
        ind = (ii-1)*M1 .+ (1:M1)
        A[:,ind] .= Diagonal(A2[:,ii]) * A1 
        ks[1:end-1,ind] .= ks1
        ks[end,ind] .= ks2[ii]
    end
    
    A, ks
end

function get_torus_err(mode_bank::AbstractArray, B_val::AbstractArray, ind_val::AbstractVector, 
                       QR::AdaptiveQR)
    kind = get_kinds(QR)
    X = get_X(QR)
    A_val = mode_bank[ind_val, kind]
    v = A_val * X - B_val
    norm(v)
end


# """
    

# """
function adaptive_get_torus(w0::AbstractVector, hs::AbstractMatrix; 
                            validation_fraction::Number=0.05, ridge_factor::Number=10)
    # Problem dimensions
    Nw0 = length(w0)
    D, N = size(hs)
    
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
    
    QRs = [AdaptiveQR(ComplexF64, N_train, D) for _ = 1:Nw0]
    for ii = 1:Nw0; 
        update!(QRs[ii], mode_bank[ind_train,ind_A],ind_A)
        solve!(QRs[ii],B_train)
    end

    dir_best = 1;
    err_best = get_torus_err(mode_bank, B_val, ind_val, QRs[dir_best]) / err_den
    M_best   = 1;

    Lmode = size(mode_bank,2);
    
    

    # Adaptively find the best validation error
    for _ = 1:N_train # This loop should always break, but making it finite for error catching
        # Consider updating in each dimension
        err_candidates = zeros(Nw0)

        for jj = 1:Nw0
            kind_update = update_inds((@view ks_bank[:,1:Lmode]), get_kinds(QRs[jj]), sz, jj)
            # println("jj=$jj, ks = ")
            # display(ks_bank[:, get_kinds(QRs[jj])])
            # println("Next ks = ")
            # display(ks_bank[:, kind_update])
            
            if length(kind_update)+QRs[jj].M > N_train
                err_candidates[jj] = Inf
            else
                update!(QRs[jj],(mode_bank[ind_train,kind_update]), kind_update)
                solve!(QRs[jj], B_train)
                err_candidates[jj] = get_torus_err(mode_bank, B_val, ind_val, QRs[jj])/err_den
            end
        end
        
        # Update in the best direction or break
        dir = argmin(err_candidates)
        err_dir = err_candidates[dir]
        # println("ks = $ks")
        # display(ks_bank)
        # println("err=$err")
        if err_dir < ridge_factor*err_best
            # If this is the best we've seen, update best values
            if err_dir < err_best
                err_best = err_dir
                dir_best = dir
                M_best   = QRs[dir].M
            end
            # println("dir = $(dir)")

            # Either way, we take a next stpe
            mode_next, ks_next = get_update_matrix(w0,N,dir,sz .+ 1)

            Lnext = size(mode_next,2)
            while Lnext+Lmode >= size(mode_bank,2)
                tmp = mode_bank
                tmp_k = ks_bank
                
                mode_bank = zeros(ComplexF64,N,2*size(tmp,2))
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
    tor = FourierTorus(D, Nks; τ = 2π .* w0);
    for (ii, k) in enumerate(eachcol(ks))
        tor.a[:, 1, (k + Nks .+ 1)...] = X[ii, :]
    end

    return tor, err_best
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
    # display(N_est)

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
    # display(length(ks))
    for (ii, k) in enumerate(eachcol(ks))
        # display(k + Nks .+ 1)
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
