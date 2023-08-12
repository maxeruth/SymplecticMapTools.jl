function periodic_f(FJ, q)
    function my_f(xs)
        F = zeros(2*q)
        xs = reshape(xs, 2, q);
        for k = 1:q
            kn = mod(k,q)+1
            Fk, Jk = FJ(xs[:,k])
            F[2*k-1:2*k] = Fk - xs[:,kn]
        end

        return F'*F/2;
    end

    return my_f
end

function periodic_g!(FJ, q)
    function my_g!(G, xs)
        xs = reshape(xs, 2, q);

        F = zeros(2q)
        J = zeros(2q, 2q)

        for k = 1:q
            kn = mod(k,q)+1
            Fk, Jk = FJ(xs[:,k])

            F[2*k-1:2*k] = Fk - xs[:,kn]
            J[2*k-1:2*k,2*k-1:2*k] = Jk
            if q == 1
                J[:,:] += Matrix(-1.0*I, 2, 2);
            else
                J[ 2*k-1:2*k , 2*kn-1:2*kn ] = Matrix(-1.0*I, 2, 2);
            end
        end
        G[:] = J'*F;
    end

    return my_g!
end

"""
    BFGS_periodic(FJ::Function, x::AbstractVector, q::Integer;
                  maxiter::Integer=50)

Find a periodic orbit using BFGS from the Optim package.

Arguments:
- `FJ`: A function that returns the map and its derivative `F, dFdx = FJ(x)`
- `x`: An initial point to find the orbit from
- `q`: The period of the orbit
- `maxiter=50`: The maximum number of optimization steps allows

Output:
- `xs`: A periodic trajectory of length `d`
- `res`: An Optim return object, (can check convergence with
  `Optim.converged(res)`)
"""
function BFGS_periodic(FJ::Function, x::AbstractVector, q::Integer;
                       maxiter::Integer=50)
    xs = zeros(d, q)
    xs[:,1] = x
    for k = d:q
        xs[:,k] = FJ(xs[:,k-1])[1]
    end
    xs = vec(xs)

    f = periodic_f(FJ, q)
    g! = periodic_g!(FJ, q)
    res = optimize(f, g!, vec(xs), BFGS(), Optim.Options(iterations=maxiter))

    return reshape(Optim.minimizer(res), d, q), res;
end

function periodic_resid!(F::AbstractArray, J::AbstractArray, FJ::Function,
                         xs::AbstractArray, q::Integer)
    # Fill F with residuals of (F(k) = x(k+1)-F(x(k)) )
    # Fill diagonal of J with map Jacobian ( J(k,k) = JF(x(k)) )
    for k = 1:q
        kn = mod(k,q)+1
        Fk, Jk = FJ(xs[:,k])
        F[d*k-1:d*k] = xs[:,kn] - Fk
        J[d*k-1:d*k,d*k-1:d*k] = Jk
        if q == 1
            J[:,:] += Matrix(-1.0*I, d, d);
        end
    end
end

"""
    newton_periodic(FJ::Function, x::AbstractVector, q::Integer;
                    maxiter::Integer=50, rtol::Number=1e-8, verbose::Bool=false)

Find a periodic orbit using Newton's method with line search

Arguments:
- `FJ`: A function that returns the map and its derivative `F, dFdx = FJ(x)`
- `x`: An initial point to find the orbit from
- `q`: The period of the orbit
- `maxiter=50`: The maximum number of optimization steps allows
- `rtol=1e-8`: The residual tolerance required
- `verbose=false`: If true, outputs information about convergence

Output:
- `xs`: A periodic trajectory of length `d`
- `converged`: A flag that indicates whether the orbit converged in `maxiter`
  steps
"""
function newton_periodic(FJ::Function, x::AbstractVector, q::Integer;
                         maxiter::Integer=50, rtol::Number=1e-8,
                         verbose::Bool=false)
    # Initialize the orbit from guess point
    xs = zeros(d, q)
    xs[:,1] = x
    for k = 2:q
        xs[:,k] = FJ(xs[:,k-1])[1]
    end

    J = zeros(d*q, d*q)
    F = zeros(d*q)
    # Initialize the Jacobian to '-I' in a block-shift-map pattern ( J(k,k+1) = -I )
    for k = 1:q
        kn = mod(k,q)+1
        J[ d*k-1:d*k , d*kn-1:d*kn ] = Matrix(-1.0*I, d, d);
    end
    periodic_resid!(F, J, FJ, xs, q)
    normF = norm(F)
    normk = 2*normF;

    for step = 1:maxiter
        dx = reshape(J\F, d, q)
        xs_prev = copy(xs);
        # Solve for Newton step, return if we are done

        α = 1.0
        for k = 1:10
           xs = xs_prev + α.*dx;
           periodic_resid!(F, J, FJ, xs, q)

           normk = norm(F);

           if (normk < normF)
              break
           end
           α /= 2
        end
        normF = normk

        if verbose
           println("newton_periodic step=$step normF=$normF α=$α")
        end

        if normF < rtol
            if verbose
                println("Converged in $step steps")
            end
            return xs, true
        end
    end

    if verbose
        println("Warning: newton_periodic did not converge in $(maxiter) steps: $(norm(F)) > $(rtol)")
    end
    return xs, false
end

function FJ_qfold(FJ, q)
    function FJq(x)
        J = [1.0 0.0 ; 0.0 1.0]
        for k = 1:q
            x, Jk = FJ(x)
            J = Jk*J
        end
        return x, J
    end

    return FJq
end
