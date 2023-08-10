include("./MPEOperators.jl")
include("./MPE.jl")
include("./ContinuedFractions.jl")
using Polynomials

function invariant_circle_model(h::Function, F::Function, x0::AbstractVector,
                                N::Integer, K::Integer; iterative::Bool=true,
                                x_prev::Union{AbstractArray,Nothing}=nothing)
    x = deepcopy(x0);
    h0 = h(x);
    d = length(h0);

    @assert N*d > 2K-1

    hs = zeros(d, N+2K-2);
    xs = zeros(2, N+2K-2);

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

    for ii = ii_init:N+2K-2
        x = F(x);
        hs[:,ii] = h(x)
        xs[:,ii] = x;
    end
    history = 0

    if iterative
        c, sums, resid, history = vector_MPE_iterative(hs, K)
    else
        c, sums, resid = vector_MPE_backslash(hs, K)
    end

    return c, reshape(sums, d, N), reshape(resid, d, N-1), xs, hs, history;
end

function adaptive_invariant_circle_model(h::Function, F::Function,
                    x0::AbstractVector; rtol::Number=1e-12, Kmax = 30,
                    Kstride=20, iterative::Bool=true, Nfactor::Integer=1)
    #
    d = length(h(x0));
    K = 20
    N = Nfactor*K ÷ d;
    # println("K=$K, N=$N")

    c, sums, resid, xs, hs, history = 0,0,Inf,nothing,nothing,0;
    rnorm = Inf;
    while (K < Kmax) && (rnorm > rtol)
        K += Kstride;
        N += (Nfactor÷d)*Kstride;
        # println("K=$K, N=$N")

        c, sums, resid, xs, hs, history = invariant_circle_TS(h, F, x0, N, K; iterative, x_prev=xs)
        rnorm = norm(resid)
    end

    return c, sums, resid, xs, hs, rnorm, K, history
end

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

function get_sum_ave(hs, c)
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
    p = zeros(Int64, N);
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

function eigenvalues_and_project(hs, c; growth_tol=1e-5)
    sum_ave =  get_sum_ave(hs, c)
    λs = roots(Polynomial(c));
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

function get_circle_coef(hs::AbstractArray, ω0::Number)
    den = denoms(ContFrac(ω0/2π));
    d, N = size(hs);

    Nmode = min(floor(Int64, N/2), maximum(den))
    Nmode = mod(Nmode, 2) == 0 ? Nmode - 1 : Nmode

    # λ = exp(2π*im * ω0);
    ωs = zeros(Nmode);
    ωs[1] = 0;
    for ii = 1:Nmode÷2
        ωs[2ii] = ii*ω0
        ωs[2ii+1] = -ii*ω0
    end

    w = WBA_weight(1, N)
    w = ones(N)
    w = w / sqrt(sum(w.^2))

    vs = Diagonal(w)*[exp(m*im*ωi) for m=0:N-1, ωi = ωs]
    rhs = Diagonal(w)*hs'

    return vs \ rhs
end


"""
    get_circle_info(hs::AbstractArray, c::AbstractArray; rattol::Number=1e-8,
                    ratcutoff::Number=1e-4, max_island_d::Integer=30)

Get a Fourier representation of an invariant circle from the observations
`hs` and the learned filter `c` (see `adaptive_invariant_circle_model` and
`invariant_circle_model` to find the filter).

Optional Arguments:
- `rattol`: Roots are judged to be rational if |ω-m/n|<rattol
- `ratcutoff`: Relative prominence needed by a linear mode to qualify as
  "important" for deciding whether the sequence is an island
- `max_island_d`: Maximum denominator considered for islands.
"""
function get_circle_info(hs::AbstractArray, c::AbstractArray;
                         rattol::Number=1e-8, ratcutoff::Number=1e-4,
                         max_island_d::Integer=30)
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
        circle_coef = get_circle_coef(hs, ω0)
    else
        # If there are islands, rerun extrapolation to find reduced filter
        d, N = size(hs);
        hs_island = reshape(hs[:, 1:(N÷Nisland)*Nisland], d*Nisland, N÷Nisland)
        zs = Vector{Any}(undef, Nisland);
        K = length(c) ÷ 2;
        Kisland = K ÷ Nisland +1;
        c_island, ~, resid, ~ = vector_MPE_iterative(hs_island, Kisland);

        λs, ωs, sum_ave, coefs = eigenvalues_and_project(hs_island, c_island)

        # Then, find coefficients of each island
        ω0 = ωs[1]*2π
        circle_coef = get_circle_coef(hs_island, ω0)
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
