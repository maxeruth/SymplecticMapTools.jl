"""
    lyapunov_exponent(x0::AbstractVector, FJ::Function, N::Integer)

Compute the lyapunov exponent of the map `FJ` after `N` iterates of the map.

Input:
- `x0`: Initial point
- `FJ`: Returns the symplectic map evaluation `F` and Jacobian `J = dF/dx`, i.e.
  has the function signature `F, J = FJ(x)`
- `N`: The number of iterations of the map to be used

Output:
- `λ`: The Lyapunov exponent `λ = 1/N * log(σ_max(dF^N(x)/dx))`
"""
function lyapunov_exponent(x0::AbstractVector, FJ::Function, N::Integer)
    x = x0;
    d = length(x);
    J = Matrix(1.0*I, d, d);
    logJnorm = 1.0;

    for ii = 1:N
        x, Jx = FJ(x);
        J = Jx*J

        Jnormx = norm(J);
        J = J ./ Jnormx;
        logJnorm = logJnorm + log(Jnormx);
    end

    σs = svdvals(J);

    return (logJnorm + log(σs[1]))/N
end



"""
    lyapunov_exponents(x0::AbstractVector, FJ::Function, N::Integer)

Compute the lyapunov exponents of the map `FJ` over `N` iterates of the map.

Input:
- `x0`: Initial point
- `FJ`: Returns the symplectic map evaluation `F` and Jacobian `J = dF/dx`, i.e.
  has the function signature `F, J = FJ(x)`
- `N`: The number of iterations of the map to be used

Output:
- `λs`: The Lyapunov exponent `λ = 1/n * log(σ_max(dF^n(x)/dx))` for `1<=n<=N`
"""
function lyapunov_exponents(x0::AbstractVector, FJ::Function, N::Integer)
    x = x0;
    d = length(x);
    J = Matrix(1.0*I, d, d);
    logJnorm = 1.0;
    λs = zeros(N);

    for ii = 1:N
        x, Jx = FJ(x);
        J = Jx*J

        Jnormx = norm(J);
        J = J ./ Jnormx;
        logJnorm = logJnorm + log(Jnormx);

        σs = svdvals(J);
        λs[ii] = (logJnorm + log(σs[1]))/ii
    end

    return λs
end
