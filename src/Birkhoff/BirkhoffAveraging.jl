include("./ExtrapolationOperators.jl")
include("./MPE.jl")
include("./ContinuedFractions.jl")

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
    xs = zeros(2, Nx);

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
- `h`: The observable function (can choose the identity as a default
  `h = (x)->x`)
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
    # println("K=$K, N=$N")

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

# This is a convolution. It could probably be done with FFTs if I really cared
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

function find_rationals(ωs::AbstractVector, λs::AbstractVector, dmax::Integer, tol::Number)
    Nisland = 1;
    ind = falses(length(ωs))
    for (ii, ω) in enumerate(ωs)
        n, d = SmallDenom(ω, tol)
        # println("ii=$ii, ω=$ω, =$ω≈$(n//d)")
        if d ≤ dmax
            Nisland = lcm(d, Nisland);
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
function palindromic_cheb_roots(c::AbstractVector)
    v = palindromic_to_chebyshev(c)
    C = chebyshev_companion_matrix(v)
    eigvals(C)
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
        λs[2ii-1:2ii] = roots(Polynomial([0.5, -μs[ii], 0.5]))
    end

    λs
end

function eigenvalues_and_project(hs, c; growth_tol=1e-5)
    sum_ave =  get_sum_ave(hs, c)
    λs = cheb_roots_to_roots(palindromic_cheb_roots(c));
    λs = λs[abs.(abs.(λs) .- 1) .< growth_tol]
    if length(λs) == 0
        @warn "No eigenvalues were found on the unit circle. This means rattol was probably too large"
    end
    d, N = size(hs);
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

function get_circle_coef(hs::AbstractArray, ω0::Number, maxNmode::Integer)
    den = denoms(ContFrac(ω0/2π));
    d, N = size(hs);

    Nmode = minimum([floor(Int, N/4), maximum(den), maxNmode])
    Nmode = mod(Nmode, 2) == 0 ? Nmode - 1 : Nmode
    N = min(N, 6*Nmode)

    # λ = exp(2π*im * ω0);
    ωs = zeros(Nmode);
    ωs[1] = 0;
    for ii = 1:Nmode÷2
        ωs[2ii] = ii*ω0
        ωs[2ii+1] = -ii*ω0
    end

    w = wba_weight(1, N)
    # w = ones(N)
    w = w / sqrt(sum(w.^2))

    vs = Diagonal(w)*[exp(m*im*ωi) for m=0:N-1, ωi = ωs]
    rhs = Diagonal(w)*hs[:,1:N]'

    return vs \ rhs
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
- `maxNmode`: Maximum number of considered for the integration

Output:
- `z`: An invariant circle of type `FourierCircle`
"""
function get_circle_info(hs::AbstractArray, c::AbstractArray;
                         rattol::Number=1e-8, ratcutoff::Number=1e-4,
                         max_island_d::Integer=30, ϵ::Number=0.0,
                         maxNmode::Integer=100)
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
    ind_rat, Nisland = find_rationals(ωs[ind], λs[ind], max_island_d, rattol)

    # Get the coefficients of the invariant circle
    ω0 = ωs[1]*2π
    circle_coef = [];
    if Nisland == 1
        # Find coefficients of circles
        circle_coef = get_circle_coef(hs, ω0, maxNmode)
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
        circle_coef = get_circle_coef(hs_island, ω0, maxNmode)
    end


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

    return z
end
