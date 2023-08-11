# Connecting orbit parameterized by
"""
    ConnectingOrbit(Na::Integer, p::Integer)

Initialize a connecting orbit consisting of `p` segments
cⱼ : [0,1] -> R². The endpoints cⱼ(0) and cⱼ(1) are both hyperbolic perodic
orbits. The connecting orbits connect these, acting as the outer boundary of an
island. The connecting orbits are parameterized by Chebyshev polynomials of
degree `Na`.

See `linear_initial_connecting_orbit` and `gn_connecting!` for functions to
initialize and compute a connecting orbit.
"""
struct ConnectingOrbit
    a::AbstractArray;    # List of chebyshev coefficients
    p::Integer;          # Period of the connecting orbit

    function ConnectingOrbit(Na::Integer, p::Integer)
        a = zeros(2*(Na+1), p)
        new(a, p);
    end
end

"""
    get_Na(c::ConnectingOrbit)

Get the polynomial degree of the connecting orbit.
"""
function get_Na(c::ConnectingOrbit)
   return (size(c.a, 1)÷2) - 1;
end

"""
    get_p(c::ConnectingOrbit)

Get the period of the connecting orbit
"""
function get_p(c::ConnectingOrbit)
    return c.p;
end

function get_a(c::ConnectingOrbit; i_p::Integer=0)
    if i_p == 0
        return c.a[:,:];
    else
        return c.a[:, i_p];
    end
end

function set_a!(c::ConnectingOrbit, a::AbstractArray; i_p::Integer=0)
    if i_p == 0
        c.a[:] = a[:];
    else
        c.a[:, i_p] = a[:];
    end
end

"""
    get_am(c::ConnectingOrbit, m::Integer; i_p::Integer=1)

Get the `m`th Chebyshev coefficient of the `i_p`th connecting orbit
"""
function get_am(c::ConnectingOrbit, m::Integer; i_p::Integer=1)
   return c.a[2m+1:2m+2, i_p];
end

"""
    set_am!(c::ConnectingOrbit, m::Integer, am::AbstractArray; i_p::Integer=1)

Set the `m`th Chebyshev coefficient of the `i_p`th connecting orbit
"""
function set_am!(c::ConnectingOrbit, m::Integer, am::AbstractArray; i_p::Integer=1)
    c.a[2m+1:2m+2, i_p] = am;
end

"""
    Base.getindex(c::ConnectingOrbit, i_A::Integer, i_p::Integer)

Wrapper for `get_am`.
"""
function Base.getindex(c::ConnectingOrbit, i_A::Integer, i_p::Integer)
    return get_am(c, i_A; i_p);
end

"""
    Base.setindex!(c::ConnectingOrbit, X::AbstractArray, i_A::Integer, i_p::Integer)

Wrapper for `set_am!`
"""
function Base.setindex!(c::ConnectingOrbit, X::AbstractArray, i_A::Integer, i_p::Integer)
    set_am!(c, i_A, X; i_p)
end

"""
    evaluate(c::ConnectingOrbit, x::AbstractArray; i_p::Integer=1)

Evaluate the `i_p`th connecting orbit at a set of points `x[j]` in [0,1]
"""
function evaluate(c::ConnectingOrbit, x::AbstractArray; i_p::Integer=1)
    Na = get_Na(c);
    Nx = length(x);

    Ts = zeros(Nx, Na+1)
    Ts[:, 1] .= 1.0
    Ts[:, 2] = vec(x);
    for k = 3:Na+1
        Ts[:, k] = 2 .* vec(x) .* Ts[:, k-1] .- Ts[:, k-2]
    end

    return reshape(c.a[:, i_p], 2, Na+1)*Ts';
end

"""
    evaluate(c::ConnectingOrbit, x::Number; i_p::Integer=1)

Evaluate the `i_p`th connecting orbit at a point `x` in [0,1]
"""
function evaluate(c::ConnectingOrbit, x::Number; i_p::Integer=1)
    return vec(evaluate(c, [x]; i_p));
end

"""
    (c::ConnectingOrbit)(x, i_p)

Wrapper for evaluate(c, x; i_p)
"""
function (c::ConnectingOrbit)(x, i_p)
    return evaluate(c, x; i_p)
end

"""
    deval(c::ConnectingOrbit, x::Number; i_p::Integer=1)

Evaluate the derivative of the `i_p`th connecting orbit at a point `x` in [0,1]
"""
function deval(c::ConnectingOrbit, x::Number; i_p::Integer=1)
    Na = get_Na(c);

    Us = zeros(3)
    Us[1] = 1;
    Us[2] = 2*x;

    deriv = Us[1] * c[1, i_p] + 2 * Us[2] * c[2, i_p];

    for i_A = 3:Na
        Us[3] = 2*x*Us[2] - Us[1];
        deriv = deriv + i_A * Us[3] * c[i_A, i_p];
        Us[1:2] = Us[2:3];
    end

    return deriv;
end

function quad_nodes(c::ConnectingOrbit, Ns::Integer)
    return -cos.(LinRange(0, π, Ns))
end

"""
    area(c::ConnectingOrbit; origin::AbstractArray=[0.0,0.0],
         i_p::Integer=1, rtol=1e-8)

Get the "area" of a connecting orbit ∫₋₁¹ (x dy/ds - y dy/ds) ds
"""
function area(c::ConnectingOrbit; origin::AbstractArray=[0.0,0.0], i_p::Integer=1, rtol=1e-8)
    J = [0.0 -1.0; 1.0 0.0];
    f = s -> deval(c, s; i_p)' * J * (c(s, i_p) - origin) / 2;
    A, err = quadgk(f, -1, 1; rtol)

    return A;
end


## Accelarated routines
# Routines that work on a grid
function get_Φ(c::ConnectingOrbit, s::AbstractVector)
    Na = get_Na(c);

    Φ = zeros(length(s), 1 + Na);
    Φ[:, 1] .= 1;
    Φ[:, 2]  = s;
    for k = 3:Na+1
        Φ[:, k] = 2 .* s .* Φ[:, k-1] .- Φ[:, k-2]
    end

    return Φ;
end

# x should be size 2×Nθ×p
function grid_eval!(x::AbstractArray, c::ConnectingOrbit, Φ::AbstractArray; i_p::Integer=1)
    reshape(x, 2, size(Φ,1))[:,:] = reshape(c.a[:, i_p], 2, get_Na(c) + 1)*Φ';
end


function connecting_parameterization(λ0, λ1)
    if (λ0 < 1) || (λ1 > 1)
        println("connecting_parameterization: λ0 < 1 or λ1 > 1, rescaling")
        nrm = sqrt(λ0^2 + λ1^2);
        λ0 = λ0/nrm
        λ1 = λ1/nrm
    end

    dλ0 = λ0 - 1;
    dλ1 = 1 - λ1;

    c2 = -(dλ0 + dλ1)*(dλ1 + dλ0*(2*dλ1-1));
    c1 = 2*(1+dλ0)*dλ1;

    f = (s) -> begin
        x = (s+1)/2;
        y = x/(dλ1 + sqrt(dλ1^2 + (dλ0^2 - dλ1^2)*x));
        2*y*(c1 + c2*y) - 1
    end

    return f
end

## Seeding functions
"""
    linear_initial_connecting_orbit(x0, x1, Na::Integer)

Initialize a ConnectingOrbit of degree `Na` using the two hyperbolic orbits `x0`
and `x1` of the form `xn = [xn0, F(xn0), ..., F^(p-1)(xn0)]`.
"""
function linear_initial_connecting_orbit(x0, x1, Na::Integer)
    p = length(x0)÷2;
    c = ConnectingOrbit(Na, p)

    for i_p = 1:p
        x0_p = x0[2*i_p-1:2*i_p];
        x1_p = x1[2*i_p-1:2*i_p];
        c[0, i_p] = (x0_p + x1_p)/2;
        c[1, i_p] = (x1_p - x0_p)/2;
    end

    return c;
end


## Newton iteration
function connecting_residv!(R::AbstractArray, J::AbstractArray,
                            c::ConnectingOrbit, FJp::Function,
                            s::AbstractVector, gs::AbstractVector,
                            Φs::AbstractArray, Φg::AbstractArray)
    Ns = length(s);
    Na = get_Na(c);

    qs = zeros(2, Ns);
    grid_eval!(qs, c, Φs)
    grid_eval!(R, c, Φg);

    id = Matrix(1.0*I, 2, 2);
    J[:, :] = kron(Φg, id);

    for i_s = 1:Ns
        F_i, dF_i = FJp(qs[:, i_s]);
        R[2i_s-1:2i_s] = R[2i_s-1:2i_s] - F_i;

        J[2i_s-1:2i_s, :] = J[2i_s-1:2i_s, :] - kron(Φs[[i_s], :], dF_i);
    end
end


function connecting_constraints(c::ConnectingOrbit, FJp::Function, ends::AbstractArray; verbose=false)
    Na = get_Na(c);
    p = get_p(c);

    id = Matrix(1.0*I, 2, 2);
    con = zeros(5, 2*(Na+1));

    # Get eigenvector information
    ~, J0 = FJp(ends[:, 1, 1]);
    eigen0 = eigen(J0);
    v0 = abs(eigen0.values[1]) > 1 ? eigen0.vectors[:, 1] : eigen0.vectors[:, 2];

    dc0 = deval(c, -1.0);
    if v0'*dc0 < 0
        v0 = -v0;
    end

    ~, J1 = FJp(ends[:, 2, 1]);
    eigen1 = eigen(J1);
    v1 = abs(eigen1.values[1]) < 1 ? eigen1.vectors[:, 1] : eigen1.vectors[:, 2];

    dc1 = deval(c, 1.0);
    if v1'*dc1 < 0
        v1 = -v1;
    end

    if verbose
        println("    v0=$(v0), v1=$(v1)")
    end

    for i_A = 1:Na+1
        ind_A = 2i_A-1:2i_A;

        # Fix the left end at x0
        con[1:2, ind_A] = (-1)^(i_A-1).*id;

        # Fix the right end at x1
        con[3:4, ind_A] = id;

        # Fix the derivatives at the two ends to be equal
        con[5, ind_A] = (i_A-1)^2 * ((-1)^(i_A) * v0 - v1)
    end

    return con
end

function bisection(g::Function, x)
    if x == -1
        return -1.0
    end

    if x == 1
        return 1.0
    end

    a = -1.;
    ga = -1.

    b = 1.;
    gb = 1.;

    while (b-a) > 1e-14
        c = (a+b)/2;
        gc = g(c);
        if gc < x
            a = c;
            ga = gc;
        else
            b = c;
            gb = gc
        end
    end

    return a;
end

function get_remaining(c::ConnectingOrbit, FJ::Function, Ns::Integer)
    println("In get_remaining")
    Na = get_Na(c);
    p = get_p(c);

    s = quad_nodes(c, Ns)
    Φ = get_Φ(c, s);
    eval_mat = kron(Φ, Matrix(1.0*I, 2, 2));

    for i_p = 2:p
        v0 = deval(c, -1; i_p=i_p-1);
        ~, J0 = FJ(c(-1, i_p-1));
        λ0 = norm(J0*v0)/norm(v0);

        v1 = deval(c, 1; i_p=i_p-1);
        ~, J1 = FJ(c(1, i_p-1));
        λ1 = norm(J1*v1)/norm(v1);

        g = connecting_parameterization(λ0, λ1);
        # println("typeof(g) = $(typeof(g))")
        println("get_remaining, i_p=$i_p, λ0=$λ0, λ1=$λ1")
        ginv = (x) -> bisection(g, x);

        Fs = zeros(2, Ns);
        for i_s = 1:Ns
            Fs[:, i_s] = FJ( c(ginv(s[i_s]), i_p-1) )[1]
        end

        set_a!(c, eval_mat\vec(Fs) ; i_p)
    end
end

function get_λ_ends(x0, x1, FJp)
    ~, J0 = FJp(x0);
    eigen0 = eigen(J0);
    λs0 = eigen0.values;
    λ0 = abs(λs0[1]) > abs(λs0[2]) ? λs0[1] : λs0[2]
println("get_λ_ends: λs0 = $λs0")
    ~, J1 = FJp(x1);
    eigen1 = eigen(J1);
    λs1 = eigen1.values;
    λ1 = abs(λs1[1]) < abs(λs1[2]) ? λs1[1] : λs1[2]

    return λ0, λ1
end

"""
    gn_connecting!(c::ConnectingOrbit, FJ::Function, ends::AbstractArray;
                   Ns = 100, maxiters = 10, rtol = 1e-8, verbose=false)

Find a connecting orbit using Gauss Newton with linesearch.

Arguments:
- `c::ConnectingOrbit`: An initial connecting orbit guess, see
  `linear_initial_connecting_orbit`
- `FJ::Function`: Function defined on R² with signature  `F(x), J(x) = FJ(x)`
  where `F` is the symplectic map and `J = dF/dx`
- `ends::AbstractArray`: A 2p × 2 matrix containing the end periodic orbits
  `[x00, x10; F(x00), F(x10); ... ; F^(p-1)(x00), F^(p-1)(x10)]`
- `Ns = 100`: Number of quadrature nodes for the least squares
- `maxiters = 10`: Maximum number of Gauss Newton steps
- `rtol = 1e-8`: Stopping tolerance
"""
function gn_connecting!(c::ConnectingOrbit, FJ::Function, ends::AbstractArray;
                        Ns = 100, maxiters = 10, rtol = 1e-8, verbose=false)
    Na = get_Na(c);
    p = get_p(c);
    FJp = FJ_qfold(FJ, p);
    Nμ = 5;

    # Get g
    λ0, λ1 = get_λ_ends(ends[:, 1], ends[:, 2], FJp);
    g = connecting_parameterization(λ0, λ1)

    # Get grid and helper matrices
    s = quad_nodes(c, Ns)
    gs = g.(s);
    Φs = get_Φ(c, s)
    Φg = get_Φ(c, gs)

    # Get residual info
    R = zeros(2Ns)
    J = zeros(2Ns, 2*(Na+1))
    rhs = zeros(2*(Na+1) + Nμ);
    H = zeros(2*(Na+1) + Nμ, 2*(Na+1) + Nμ);

    if verbose
println("Starting gn_connecting! iteration with Ns = $(Ns), p = $(p)")
println("    unstable_end = $(ends[:, 1]), stable_end = $(ends[:, 2])")
    end

    # Get constraint
    aind = 1:2*(Na+1)
    μind = 2*(Na+1) .+ (1:Nμ);
    con = connecting_constraints(c, FJp, ends; verbose)

    connecting_residv!(R, J, c, FJp, s, gs, Φs, Φg)
    rhs[aind] = -J'*R;
    rhs[μind[1:2]] = ends[:, 1] - con[1:2, :] * get_a(c; i_p=1);
    rhs[μind[3:4]] = ends[:, 2] - con[3:4, :] * get_a(c; i_p=1);
    rhs[μind[5]] = - con[5, :]' * get_a(c; i_p=1);

    rnorm = norm(rhs)/Ns

    for k = 1:maxiters
        H[aind, aind] = J'*J;
        H[μind, aind] = con
        H[aind, μind] = con'

        # rhs[aind] = -J'*R;
        # rhs[μind[1:2]] = ends[:, 1] - con[1:2, :] * get_a(c; i_p=1);
        # rhs[μind[3:4]] = ends[:, 2] - con[3:4, :] * get_a(c; i_p=1);

        da = H\rhs;
        a = vec(get_a(c; i_p=1))

        α = 1.0
        for k = 1:10
            set_a!(c, a + α.*da[aind]; i_p=1)
            connecting_residv!(R, J, c, FJp, s, gs, Φs, Φg)

            rhs[aind] = -J'*R;
            rhs[μind[1:2]] = ends[:, 1] - con[1:2, :] * get_a(c; i_p=1);
            rhs[μind[3:4]] = ends[:, 2] - con[3:4, :] * get_a(c; i_p=1);
            rhs[μind[5]]   = - con[5, :]' * get_a(c; i_p=1);

            normk = norm(rhs)/Ns
            if (normk < rnorm) || (k == 10)
                rnorm = normk;
                break
            end
            α /= 2
        end

        rnorm = norm(rhs)/Ns
        if verbose
println("    k=$k, rnorm=$rnorm, α=$α")
        end
        if rnorm < rtol
            get_remaining(c, FJ, Ns)
            return rhs
        end
    end

    println("Did not converge in $(maxiters) steps: $(rnorm) > $(rtol)")
    get_remaining(c, FJ, Ns)
    return rhs
end
