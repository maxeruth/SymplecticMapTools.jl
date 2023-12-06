function get_sums_and_resid(x,c,D,H)
    N = length(c);
    X = block_hankel_linear_operator(x[:, 1:end-N+1], x[:, end-N+1:end])
    sums = X*c;
    resid = D*H*c;

    sums, resid
end

"""
    vector_mpe_backslash(x::AbstractArray, K::Integer)

Applies Birkhoff vector MPE to a sequence `x_n = x[:, n]`
"""
function vector_mpe_backslash(x::AbstractArray, K::Integer; ϵ = 0.0)
    if ϵ != 0.0; println("Warning: vector_mpe_iterative not supported for ϵ!=0.0"); end
    x = typeof(x) <: AbstractVector ? x' : x
    u = diff(x, dims=2)
    d, L = size(u);

    N = 2K-1;
    M = L - N - 1;
    @assert (d*M ≥ K)

    P = Matrix(mpe_p(K-1));
    D = wba_weight(d, M);
    H = BlockHankelMatrix(u[:,2:M+1], u[:,M+1:N+M])

    A = D*H*P';
    b = - D*(vec(u[:,1:M]) + vec(u[:,end-M+1:end]));

    c = ones(2K+1);
    c[2:end-1] = (P'*(A\b))
    c = c ./ sum(c)

    H2 = BlockHankelMatrix(u[:,1:end-2K], u[:,end-2K:end])
    sums, resid = get_sums_and_resid(x,c,D,H2)

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
                              btol = 1e-14, c0 = nothing, ϵ = 0.0)
    if ϵ != 0.0; println("Warning: vector_mpe_iterative not supported for ϵ!=0.0"); end
    x = typeof(x) <: AbstractVector ? x' : x
    u = diff(x, dims=2)
    d, L = size(u);

    N = 2K+1;
    M = L - N + 1;
    @assert (d * M ≥ K)


    D = LinearOperator(wba_weight(d, M); symmetric=true, hermitian=true);
    H = block_hankel_linear_operator(u[:,1:end-N+1], u[:,end-N+1:end]);
    P = mpe_p_2(K);

    A = D*H*P';
    b = - D*H[:, K+1]*[1.];

    c = (c0==nothing) ? ones(2K+1) : c0;
    c = c ./ c[K+1]

    sol = c[K+2:end];
    sol, history = lsqr!(sol, A, b; atol, btol, log=true)
    c[K+2:end] = sol
    c[1:K] = sol[end:-1:1];
    c = c ./ sum(c)

    sums, resid = get_sums_and_resid(x,c,D,H)

    return c, sums ,resid, history
end

# Iterative MPE with dense matrices. Used for testing
function vector_mpe_iterative_full(x::AbstractArray, K::Integer; atol = 1e-14, btol = 1e-14, ϵ=0.0)
    if ϵ != 0.0; println("Warning: vector_mpe_iterative_full not supported for ϵ!=0.0"); end
    x = typeof(x) <: AbstractVector ? x' : x
    u = diff(x, dims=2)
    d, L = size(u);

    N = 2K-1;
    M = L - N - 1;
    # @assert (d*M ≥ K)

    D = wba_weight(d, M);
    H = BlockHankelMatrix(u[:,2:M+1], u[:,M+1:N+M]);
    P = Matrix(mpe_p(K-1));


    display(size(D))
    display(size(H))
    display(size(P'))
    A = D*H*P';
    b = - D*(vec(u[:,1:M]) + vec(u[:,end-M+1:end]));

    c = ones(2K+1);
    sol, history = lsqr(A,b; atol, btol, log=true)
    c[2:end-1] = P'*sol
    c = c ./ sum(c)

    H2 = BlockHankelMatrix(u[:,1:end-2K], u[:,end-2K:end])
    sums, resid = get_sums_and_resid(x,c,D,H2)

    return c, sums ,resid, history
end

function vector_rre_iterative(x::AbstractArray, K::Integer; atol=1e-14, btol=1e-14, ϵ=0.0)
    if ϵ != 0.0; println("Warning: vector_rre_iterative not supported for ϵ!=0.0"); end
    x = typeof(x) <: AbstractVector ? x' : x
    u = diff(x, dims=2)
    d, L = size(u);

    N = 2K+1;
    M = L - N + 1;
    @assert (d * M ≥ K)

    D = LinearOperator(wba_weight(d, M); symmetric=true, hermitian=true);
    H = block_hankel_linear_operator(u[:,1:end-N+1], u[:,end-N+1:end]);
    P = rre_p(K);

    A = D*H*P';
    b = - D*H[:, K+1]*[1.];

    xi = ones(K)
    # println("MPE: size(A) = $(size(A))")
    sol, history = lsqr!(xi, A, b; atol, btol, log=true)
    c = P'*xi;
    c[K+1] += 1;

    sums, resid = get_sums_and_resid(x,c,D,H)

    return c, sums, resid, history
end


function vector_rre_backslash(x::AbstractArray, K::Integer; atol=1e-14, btol=1e-14, ϵ=0.0)
    x = typeof(x) <: AbstractVector ? x' : x
    u = diff(x, dims=2)
    d, L = size(u);

    N = 2K+1;
    M = L - N + 1;
    @assert (d * M ≥ K)

    P = Matrix(rre_p(K));
    D = wba_weight(d, M);
    H = BlockHankelMatrix(u[:,1:end-N+1], u[:,end-N+1:end]);

    A = D*H*P';
    b = - D*H[:, K+1];
    if ϵ != zero(typeof(ϵ))
        WKinv = sqrt(ϵ)*inv(wba_weight(1, 2K+1))
        A = vcat(A, WKinv*P')
        b = vcat(b, WKinv[:,K+1])
    end

    c = P'*(A\b);
    c[K+1] += 1;

    sums, resid = get_sums_and_resid(x,c,D,H)

    return c, sums, resid
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
