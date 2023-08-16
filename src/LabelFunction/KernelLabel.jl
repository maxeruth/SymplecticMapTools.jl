include("./Kernels/kernels.jl")

function constraint_matrix(N)
    return sparse(1:2N, kron(1:N, [1,1]), kron(ones(N), [1.,-1.]))
end

"""
    kernel_sample_F(F::Function, N::Integer, xb::AbstractVector,
                    yb::AbstractVector)

Sobol sample `N` points in the rectangle `xb` × `yb`. Then, evaluate `F` at
each point. Input can be used for `kernel_eigs` and `kernel_bvp`

Arguments:
- `F`: The symplectic map from the state space to itself
- `N`: The number of points to be sampled
- `xb` × `yb`: The 2-vector bounds of the rectangle to be sampled

Output:
- `xs`: A `2` × `2N` array of samples compatable with `kernel_eigs` and
  `kernel_bvp`, of the form `[x_1, F(x_1), x_2, F(x_2), ...]`
"""
function kernel_sample_F(F::Function, N::Integer, xb::AbstractVector,
                         yb::AbstractVector)
    xs = zeros(2, 2, N);
    s = SobolSeq([xb[1], yb[1]], [xb[2], yb[2]])

    for ii = 1:N
        x = next!(s)
        Fx = F(x);

        xs[:, 1, ii] = x;
        xs[:, 2, ii] = Fx;
    end

    return reshape(xs[:,:,:], 2, 2N);
end

"""
    window_weight(xs::AbstractArray, lims::AbstractVector, α::Number;
                       ind::Integer=2)

A simple boundary weighting function for `kernel_eigs` and `kernel_bvp`. Returns
a weight that is approximately 1 outside of `lims` and 0 inside via the sum of
two logistic functions.

Arguments:
- `xs`: A `d` × `N` array of points, where `d` is the size of the phase space
  and `N` is the number of points.
- `lims`: A 2-vector giving the interval where the weight is approximately 0
- `α`: The length scale of the logistic functions (typically small relative to
  the size of the domain)
- `ind=2`: The index over which the window is applied. Defaults to `2` for
  maps on T×R (such as the standard map)

Output:
- `w`: A `N`-vector with the window weights
"""
function window_weight(xs::AbstractArray, xlims::AbstractVector, α::Number;
                       ind::Integer=2)
    f = (x) -> (1. / (1. + exp((xlims[2] - x[ind])/α)) +
                1. / (1. + exp(-(xlims[1] - x[ind])/α)))
    [f(x) for x = eachcol(xs)]
end

"""
    rectangular_window_weight(xs::AbstractArray, xlims::AbstractVector,
                              ylims::AbstractVector, α::Number)

A simple boundary weighting function for `kernel_eigs` and `kernel_bvp`. Returns
a weight that is approximately 1 outside of `xlims` × `ylims` and 0 inside via
the a function of logistic functions.

Arguments:
- `xs`: A `d` × `N` array of points, where `d` is the size of the phase space
  and `N` is the number of points.
- `xlims` and `ylims`: 2-vectors giving the rectangle where the weight is
  approximately 0
- `α`: The length scale of the logistic functions (typically small relative to
  the size of the domain)

Output:
- `w`: A `N`-vector with the window weights
"""
function rectangular_window_weight(xs::AbstractArray, xlims::AbstractVector,
                                   ylims::AbstractVector, α::Number)
    w1 = window_weight(xs, xlims, α; ind=1)
    w2 = window_weight(xs, ylims, α; ind=2)
    return 1. .- (1. .- w1).*(1. .- w2)
end


"""
    kernel_eigs(xs::AbstractArray, ϵ::Number, nev::Integer, σ::Number,
                boundary_weights::Vector; kernel::Symbol=:SquaredExponential,
                zero_mean = false, check = 1)

Solve the the invariant eigenvalue problem, given by the Rayleigh quotient\\
> `min_c (‖GKc‖² + ‖Kc‖²_bd + ϵ‖c‖²_k)/‖Kc‖²` \\
where
- `‖GKc‖²` is a norm penalizing invariance (and possible a non-zero mean)
- `‖Kc‖²_bd` is a norm penalizing boundary violation
- `‖c‖²_k` is the smoothing kernel norm
- `‖Kc‖²` is the ℓ² norm of the points
The eigenvalue problem is solved via `Arpack.jl`

Arguments:
- `xs`: interpolation points of size d × 2N, where xs[:, N+1:2N] = F.(xs[:, 1:N])
- `ϵ`: Amount of regularization
- `nev`: Number of eigenvalues to find
- `σ`: Kernel width
- `boundary_weights`: Boundary weighting vector, should be positive and O(1) at
  points `x` where one wants |k(x)| << 1
- `kernel`: Type of kernel to interpolate (see `KernelLabel`)
- `zero_mean = false`: Set to true to add a constraint that encourages `k` to
  have a zero mean. Useful when `xs` are sampled on an invariant set and
  boundary_weights=0
- `check = 1`: See `Arpack.eigs`. If 1, return all of the converged eigenvectors
  and eigenvalues. If 0, throw an error instead.

Output:
- `λs`: The eigenvalues
- `vs`: The eigenvectors
- `k`: The kernel label. Use `set_c!(k, vs[:,n])` to load the `n`th eigenvector
  into `k`. By default, `k` stores the lowest order eigenvector.
"""
function kernel_eigs(xs::AbstractArray, ϵ::Number, nev::Integer, σ::Number,
                     boundary_weights::Vector; kernel::Symbol=:SquaredExponential,
                     zero_mean = false, check = 1)
    d = size(xs, 1);
    N = size(xs, 2)÷2;

    if (size(xs, 2) - 2N != 0) # Need to use an even number points
        xs = xs[:, 1:2N]
    end

    k = KernelLabel(xs, zeros(2N), σ; kernel)

    # Obtain and factorize kernel matrix
    K = Symmetric(get_matrix(k, xs));
    chol = cholesky(K, RowMaximum(), check=false)
    p = chol.piv;
    ip = invperm(p);
    L = chol.L[:, 1:chol.rank];

    # Get constraint matrix
    i1 = 1:2:2N; i2 = 2:2:2N;
    GPtL_inv = L[ip[i1], :] - L[ip[i2], :]
    v_mean = ones(2N) ./ sqrt(2N);
    GPtL_mean = zero_mean ? (ones(2N)'*L ./ sqrt(2N))' : zeros(0, chol.rank)

    # Find Eigenvalue problem matrices
    A = Symmetric(L'*L)

    GPtL = vcat(GPtL_inv, GPtL_mean);
    W = Diagonal(boundary_weights)
    B = Symmetric(GPtL'*GPtL + ϵ*I + L[ip, :]'*W*L[ip, :])

    # Call eigs
    λs, vs = Arpack.eigs(A,B; nev = nev, which = :LR, check = check);

    # Post process
    λs = 1 ./ λs;

    vs = L*(A \ vs) # Can we avoid this?
    vs = vs[ip, :]
    vs = vs * Diagonal([1. ./ norm(K*v) for v = eachcol(vs)])

    Nλ = length(λs)
    if Nλ ≥ 1
        set_c!(k, vs[:, 1])
    end

    if Nλ == nev # Arpack converged
        return λs, vs, k
    end

    # Arpack did not converge
    λs_new = zeros(nev);
    vs_new = zeros(size(vs,1), nev);
    λs_new[1:Nλ] = λs;
    λs_new[Nλ+1:end] .= NaN;
    vs_new[:, 1:Nλ] = vs;
    vs_new[:, Nλ+1:end] .= NaN;

    return λs_new, vs_new, k
end


"""
    kernel_bvp(xs::AbstractArray, ϵ::Number, σ::Number,
               boundary_weights::AbstractVector,
               boundary_values::AbstractVector;
               kernel::Symbol=:SquaredExponential, residuals::Bool=true)

Solve the the invariant boundary value least-squares problem\\
> `min_c ‖GKc‖² + ‖Kc - h_{bd}‖²_bd + ϵ‖c‖²_k = R_inv + R_bd + R_eps`\\
where
- `‖GKc‖²` is a norm penalizing invariance (and possible a non-zero mean)
- `‖Kc - h_{bd}‖²_bd` is a norm penalizing the function from violating the
  boundary condition
- `‖c‖²_k` is the smoothing kernel norm

Arguments:
- `xs`: interpolation points of size d × 2N, where xs[:, N+1:2N] = F.(xs[:, 1:N])
- `ϵ`: Amount of regularization
- `σ`: Kernel width
- `boundary_weights`: Length `2N` boundary weighting vector, should be positive
  and O(1) at points `x` where one wants |k(x)| << 1
- `boundary_values`: Length `2N` boundary value vector, indicating the value the
  function should take at each point
- `kernel=:SquaredExponential`: Type of kernel to interpolate
  (see `KernelLabel`)
- `residuals=true`: True if you want the problem to return residuals.

Output:
- `k`: The kernel function
(if `residuals=true`)
- `R`: The total residual
- `R_bd`: The boundary condition residual
- `R_inv`: The invariance residual
- `R_eps`: The smoothness residual
"""
function kernel_bvp(xs::AbstractArray, ϵ::Number, σ::Number,
                    boundary_weights::AbstractVector,
                    boundary_values::AbstractVector;
                    kernel::Symbol=:SquaredExponential, residuals::Bool=true)
    d = size(xs, 1);
    N = size(xs, 2)÷2;

    if (size(xs, 2) - 2N != 0) # Need to use an even number
        xs = xs[:, 1:2N]
    end

    if σ == 0.
        xa = [minimum(xs[ii, :]) for ii = 1:d]
        xb = [maximum(xs[ii, :]) for ii = 1:d]
        σ  = 5. * (prod(xb-xa)/(2N))^(1/d) # Arbitrary choice assuming evenly distributed points
    end

    k = KernelLabel(xs, zeros(2N), σ; kernel)

    K = Symmetric(get_matrix(k, xs));

    A = zeros(2N, 2N);
    i1 = 1:2:2N; i2 = 2:2:2N;

    # Enforce invariance constraint
    GK = K[i1, :] - K[i2, :];
    A[i1, :] = GK
    A[i2, :] = -GK;

    # Boundary and regularization
    W = Diagonal(boundary_weights)
    A = A + W * K + ϵ*I;

    # rhs
    rhs = W * boundary_values;

    # Solve
    set_c!(k, A\rhs);


    if residuals
        G = constraint_matrix(N)
        c = get_c(k);
        tmp = K*c;
        tmp = tmp[i1]-tmp[i2]
        R_inv = tmp'*tmp;

        tmp = K*c-boundary_values;
        R_bd = tmp'*W*tmp;

        R_eps = ϵ*(c'*K*c);

        R = R_bd + R_inv + R_eps;

        return k, R, R_bd, R_inv, R_eps
    else
        return k
    end
end


"""
    get_energies(k::KernelLabel; W = 0. * I)

Get the relevant energies of a kernel label `k` with boundary weighting
matrix `W`.

Output energies:
- `EK = ‖c‖²_k` is the smoothing kernel norm
- `EInv = ‖GKc‖²` is a norm penalizing invariance
- `Ebd = ‖Kc‖²_W` is a norm penalizing boundary violation
- `EL2 = ‖Kc‖²` is the ℓ² norm of the points
"""
function get_energies(k::KernelLabel; W = 0. * I)
    N = get_N(k)
    K = get_matrix(k, get_x(k));
    G = constraint_matrix(N)
    c = get_c(k)

    Kc = K*c;
    GKc = Kc[1:2:end] - Kc[2:2:end];

    EK = c'*Kc;
    EInv = GKc'*GKc;
    Ebd = Kc'*W*Kc;
    EL2 = Kc'*Kc;

    return EK, EInv, Ebd, EL2
end



"""
    kernel_birkhoff(xs::AbstractArray, fs::AbstractVector, ϵ::Number, μ::Number;
                    σ::Number=0., kernel::Symbol=:SquaredExponential,
                    boundary_points::AbstractVector = [])

This function is mostly deprecated. The solutions to the infinitely discretized
limit of this problem (the Birkhoff average) don't live in the native space. So,
the results can be odd, hard to interpret, and wiggly. Use at your own risk.

Find the "Birkhoff average" of a function using the kernel approach. Solves the
least-squares problem\\
> `min_c μ⁻¹(‖GKc‖² + ‖Kc‖²_bd) + ϵ‖c‖²_k + ‖Kc - f‖²`\\
where
- `‖GKc‖²` is a norm penalizing invariance
- `‖c‖²_k` is the smoothing kernel norm
- `‖Kc - f‖²` is a least-squares fit norm
- `‖Kc‖²_bd` is a norm penalizing boundary violation

Arguments:
- `xs`: interpolation points of size d × 2N, where xs[:, N+1:2N] = F.(xs[:, 1:N])
- `fs`: function values at points of size N
- `ϵ`: Amount of regularization (can be set to zero)
- `μ`: Weighting of invariance penalization to fit (should be small, e.g. 1e-6)
- `σ`: Scale of the problem
- `kernel`: Type of kernel to interpolate (see `KernelLabel`)
- `boundary_points`: A list of indices of points on the boundary

Output:
- `k`: A KernelLabel object
"""
function kernel_birkhoff(xs::AbstractArray, fs::AbstractVector, ϵ::Number,
                         μ::Number, σ::Number; kernel::Symbol=:SquaredExponential,
                         boundary_points::AbstractVector = [])
    d = size(xs, 1);
    N = size(xs, 2)÷2;

    if (size(xs, 2) - 2N != 0) # Need to use an even number
        xs = xs[:, 1:2N]
    end

    k = KernelLabel(xs, zeros(2N), σ; kernel)

    K = Symmetric(get_matrix(k, xs));
    G = constraint_matrix(N)

    # LHS Matrix
    A = zeros(2N, 2N);

    # Enforce invariance constraint
    GK = K[1:2:2N, :] - K[2:2:2N, :];
    A[1:2:2N, :] = GK
    A[2:2:2N, :] = -GK;

    # Enforce zero Dirichlet BC
    bd = boundary_points
    A[boundary_points, :] = A[boundary_points, :] + K[boundary_points, :];

    # Add interpolation and regularization terms
    A = (1/μ)*A + K + ϵ*I

    # Solve
    set_c!(k, A\fs);

    return k
end
