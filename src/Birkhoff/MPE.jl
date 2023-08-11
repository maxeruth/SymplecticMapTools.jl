"""
    vector_mpe_backslash(x::AbstractArray, K::Integer)

Applies Birkhoff vector MPE to a sequence `x_n = x[:, n]`
"""
function vector_mpe_backslash(x::AbstractArray, K::Integer)
    x = typeof(x) <: AbstractVector ? x' : x
    u = diff(x, dims=2)
    d, L = size(u);

    N = 2K-3;
    M = L - N - 1;
    @assert (d*M ≥ K)

    P = Matrix(mpe_p(K-1));
    D = wba_weight(d, M);

    A = D*BlockHankelMatrix(u[:,2:M+1], u[:,M+1:N+M])*P';
    b = - D*(vec(u[:,1:M]) + vec(u[:,end-M+1:end]));

    c = ones(2K-1);
    c[2:end-1] = (P'*(A\b))
    c = c ./ sum(c)

    X = BlockHankelMatrix(x[:, 1:end-2K+2], x[:, end-2K+2:end])
    sums = X*c;
    resid = A*c[K:end-1] - b*c[1];

    return c, sums, resid
end

"""
    vector_mpe_backslash(x::AbstractArray, K::Integer)

Applies Birkhoff vector MPE to a sequence `x_n = x[:, n]` using the LSQR
algorithm. This currently does not have preconditioning, and therefore is less
accurate than `vector_mpe_backslash`.

Arguments:
- `x`: The sequence
- `K`: The number of unknowns in the filter
- `c0`: The initial guess of
- `atol`, `btol`: Tolerances. See `IterativeSolvers.lsqr!`
"""
function vector_mpe_iterative(x::AbstractArray, K::Integer; atol = 1e-14,
                              btol = 1e-14, c0 = nothing)
    x = typeof(x) <: AbstractVector ? x' : x
    u = diff(x, dims=2)
    d, L = size(u);

    N = 2K-1;
    M = L - N + 1;
    @assert (d * M ≥ K)

    P = mpe_p_2(K);
    D = LinearOperator(wba_weight(d, M); symmetric=true, hermitian=true);

    H = block_hankel_linear_operator(u[:,1:end-N+1], u[:,end-N+1:end]);
    A = D*H*P';
    b = - D*H[:, K]*[1.];

    c = (c0==nothing) ? ones(2K-1) : c0;
    c = c ./ c[K]

    sol = c[K+1:end];
    # println("MPE: size(A) = $(size(A))")
    sol, history = lsqr!(sol, A, b; atol, btol, log=true)
    c[K+1:end] = sol
    c[1:K-1] = sol[end:-1:1];
    c = c ./ sum(c)

    X = block_hankel_linear_operator(x[:, 1:end-2K+2], x[:, end-2K+2:end])
    sums = X*c;
    resid = D*H*c;

    return c, sums ,resid, history
end

# Iterative MPE with dense matrices. Used for testing
function vector_mpe_iterative_full(x::AbstractArray, K::Integer; atol = 1e-14, btol = 1e-14)
    x = typeof(x) <: AbstractVector ? x' : x
    u = diff(x, dims=2)
    d, L = size(u);

    N = 2K-3;
    M = L - N - 1;
    # @assert (d*M ≥ K)

    P = Matrix(mpe_p(K-1));
    D = wba_weight(d, M);

    A = D*BlockHankelMatrix(u[:,2:M+1], u[:,M+1:N+M])*P';
    b = - D*(vec(u[:,1:M]) + vec(u[:,end-M+1:end]));

    c = ones(2K-1);
    sol, history = lsqr(A,b; atol, btol, log=true)
    c[2:end-1] = P'*sol
    c = c ./ sum(c)

    X = BlockHankelMatrix(x[:, 1:end-2K+2], x[:, end-2K+2:end])
    sums = X*c;
    resid = A*c[K:end-1] - b*c[1];

    return c, sums ,resid, history
end

## TODO
# function vector_MPE_preconditioned(x::AbstractArray, K::Integer; atol = 1e-14, btol = 1e-14, c0 = nothing)
#     x = typeof(x) <: AbstractVector ? x' : x
#     u = diff(x, dims=2)
#     d, L = size(u);
#
#     N = 2K-1;
#     M = L - N + 1;
#     @assert (d*M ≥ K)
#
#     P = mpe_p_2(K);
#     D = LinearOperator(wba_weight(d, M); symmetric=true, hermitian=true);
#     D2 = LinearOperator(wba_weight(1, 2K-1); symmetric=true, hermitian=true);
#     D2inv = LinearOperator(wba_weight(1, 2K-1); symmetric=true, hermitian=true);
#
#     H = block_hankel_linear_operator(u[:,1:end-N+1], u[:,end-N+1:end]);
#     A = D*H*D2*P';
#     b = - D*H[:, K]*D2[K,K]*[1.];
#
#     c = (c0==nothing) ? D2*ones(2K-1) : D2inv*c0;
#     c = c ./ c[K]
#
#     sol = c[K+1:end];
#     sol, history = lsqr!(sol, A, b; atol, btol, log=true)
#     c[K+1:end] = sol
#     c[1:K-1] = sol[end:-1:1];
#     c = D2*c;
#     c = c ./ sum(c)
#
#     X = block_hankel_linear_operator(x[:, 1:end-2K+2], x[:, end-2K+2:end])
#     sums = X*c;
#     # resid = A*c[K:end-1] - b*c[1];
#     resid = D*H*c;
#
#     return c, sums, resid, history
# end
