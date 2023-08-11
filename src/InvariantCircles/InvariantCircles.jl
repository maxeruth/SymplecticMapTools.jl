"""
    abstract type InvariantCircle

An abstract type representing a chain of invariant circles zᵢ where zᵢ₊₁ = F(zᵢ)
for i<=p and z₀ = F(zₚ₋₁). The main implementation of this type of FourierCircle,
but we are leaving it open for the potential implementations of circles.

Implementations should have the following functions for immediate portability:

Basic stuff
- Initializer
    `InvariantCircleImplementation(stuff)`
- Similar
    `Base.similar(z::InvariantCircle)`

Getters and Setters\\
 - Get number of unknown variables per circle
    `get_N(z::InvariantCircle)`\\
 - Get period of circle
    `get_p(z::InvariantCircle)`
 - Get parameters (except τ)
    `get_a(z::InvariantCircle)`
 - Set parameters (except τ)
    `set_a!(z::InvariantCircle, a::AbstractArray)`
 - Get rotation number τ in [0,2π)
    `get_τ(z::InvariantCircle)`
 - Set τ
    `set_τ!(z::InvariantCircle, τ::Number)`

Evaluation related routines
- Evaluate the circle (see `(z::InvariantCircle)(θ, i_circle)`)
    `evaluate(z::InvariantCircle, θ::AbstractVector; i_circle::Integer=1)`
- Evaluate the derivative of the circle w.r.t. θ
    `deval(z::InvariantCircle, θ::AbstractVector; i_circle::Integer=1)`
- Get a basis Φ for z evaluated at Nθ equispaced points (requires linearity)
    `get_Φ(z::InvariantCircle, Nθ::Integer; τ::Number=0.0)`
- Evaluate the invariant circle on Φ
    `function grid_eval!(x::AbstractArray, z::InvariantCircle, Φ::AbstractArray;`
                        `i_circle=0)`
- Evaluate the derivative of the invariant circle on Φ
    `grid_deval!(dx::AbstractArray, z::InvariantCircle, Φ::AbstractArray;`
                `i_circle=0)`

Constraints (for Newton iteration)
- Constrain the value at 0
    `fixed_θ0_constraint(z::InvariantCircle)`

Other routines (useful for continuation)
- Get average radius
    `LinearAlgebra.average_radius(z::InvariantCircle)`
- Get the area of the invariant circle
    `area(z::InvariantCircle; i_circle=1, Ns = 100)`
"""
abstract type InvariantCircle end
include("./FourierCircle.jl")

"""
    evaluate(z::InvariantCircle, θ::Number; i_circle::Integer=1)

Evaluate the `i_circle`th circle in `z` at the point `θ`
"""
function evaluate(z::InvariantCircle, θ::Number; i_circle::Integer=1)
   return evaluate(z,[θ]; i_circle);
end

"""
    deval(z::InvariantCircle, θ::Number; i_circle::Integer=1)

Evaluate the derivative of the `i_circle`th circle in `z` at the point `θ`
"""
function deval(z::InvariantCircle, θ::Number; i_circle::Integer=1)
   return deval(z,[θ]; i_circle, nderiv)
end

"""
    shifted_eval(z::InvariantCircle, θ::AbstractVector; i_circle::Integer=1)

Evaluate the `i_circle`th circle in `z` at the point `θ`, shifted by `τ`
"""
function shifted_eval(z::InvariantCircle, θ::AbstractVector; i_circle::Integer=1)
   return evaluate(z, θ+get_τ(z); i_circle);
end

"""
    (z::InvariantCircle)(θ, i_circle)

Wrapper for `evaluate(z::InvariantCircle, θ::Number; i_circle::Integer=1)`
"""
function (z::InvariantCircle)(θ, i_circle)
    return evaluate(z, θ; i_circle)
end

# This seems like overkill
function fixed_τ_constraint(z::InvariantCircle)
   N = get_N(z);
   c = zeros(1, N+1);
   c[1, 1] = 1.0;

   return c
end



## Residual calculations
# Residual is size R = 2×Nθ×p
# J10 is size J10 = 2×Nθ(×1)
# J11, Jii are size 2×Nθ×4Na+2
# J_off_diag is size 2×Nθ×4Na+2×p
function resid!(R, J10, J11, Jii, J_off_diag, FJ, z, Φ, Φτ)
   # Get parameters
   Nθ = size(Φ, 1);
   p  = get_p(z);
   N = get_N(z);
   id = Matrix(I*1.0, 2, 2);

   # Initialize
   J10[:,:] .= 0.0;

   x = zeros(2, Nθ, p);
   grid_eval!(x, z, Φ);

   # Fill residual with evaluations of z
   x_τ = zeros(2, Nθ)
   grid_eval!(x_τ, z, Φτ; i_circle=1);
   R[:, :, 1] = x_τ;
   if p>1
      R[:, :, 2:end] = x[:, :, 2:end]
   end

   # variation in τ
   grid_deval!(J10, z, Φτ; i_circle=1)

   # variation in δC₁(θ+τ)
   J11[:, :, :] = reshape(kron(Φτ, id), 2, Nθ, N); # It's weird this works out

   # Variation in C (only for p > 1, could be precomputed)
   if p>1
      Jii[:, :, :] = reshape(kron(Φ, id), 2, Nθ, N);
   end

   # variation of F(C(i-1)) and residual
   for i_p = 1:p
      for i_θ = 1:Nθ
         F, dF = FJ(x[:, i_θ, i_p])

         # Adjust residual
         i_R = mod1(i_p+1, p);
         R[:, i_θ, i_R] += -F;

         if p == 1
            # Adjust J11
            J11[:, i_θ, :] += kron(Φ[[i_θ], :], -dF)
         else
            # Adjust J_off_diag
            J_off_diag[:, i_θ, :, i_p] = kron(Φ[[i_θ], :], -dF)
         end
      end
   end

   return R, J10, J11, Jii, J_off_diag
end

# This algorithm fills a dense H. We can get fancy with sparsity and solvers later
function hessian!(H, J10, J11, Jii, J_off_diag, z)
   p = get_p(z);
   N = get_N(z);
   Nθ = size(J11, 2);
   quad_weight = 2π/Nθ;

   J = full_J(J10, J11, Jii, J_off_diag, z)
   H[1:end-2,1:end-2] = quad_weight.*J'*J;

   # H[1, 1] = vec(J10)'*vec(J10) * quad_weight;
   # H[2:N+1, 1] = reshape(J11, 2Nθ, N)'*vec(J10) .* quad_weight
   # H[1, 2:N+1] = H[2:N+1, 1];
   #
   # if p == 1
   #    H[2:N+1, 2:N+1] = reshape(J11, 2Nθ, N)' * reshape(J11, 2Nθ, N) .* quad_weight;
   # elseif p == 2
   #    ind0 = [1];
   #    ind1 = 2:N+1;
   #    ind2 = N+2:2N+1;
   #
   #    J10 = reshape(J10, 2Nθ, 1);
   #    J11 = reshape(J11, 2Nθ, N);
   #    J22 = reshape(Jii, 2Nθ, N)
   #    J21 = reshape(J_off_diag[:,:,:,1], 2Nθ, N);
   #    J12 = reshape(J_off_diag[:,:,:,2], 2Nθ, N);
   #
   #    H[ind2, ind0] = (J12' * J10) .* quad_weight;
   #    H[ind0, ind2] = H[ind2, 1]';
   #    H[ind1, ind1] = (J11'*J11 + J21'*J21) .* quad_weight;
   #    H[ind2, ind2] = (J22'*J22 + J12'*J12) .* quad_weight;
   #    H[ind2, ind1] = (J12'*J11 + J22'*J21) .* quad_weight;
   #    H[ind1, ind2] = H[ind2, ind1]';
   # else
   #    # Other J10 term
   #    H[end-N+1:end, 1] = vec(J10)' * reshape(J_off_diag[:,:,:,end], 2Nθ, N) .* quad_weight;
   #    H[1, end-N+1:end] = H[end-N+1:end, 1]
   #    # Diagonal blocks
   #    for i_p = 1:p
   #       ind_i_p = 2+(i_p-1)*N:1+i_p*N
   #       J_mid = reshape((i_p==1) ? J11 : Jii, 2Nθ, N);
   #       J_next = reshape(J_off_diag[:,:,:,i_p], 2Nθ, N);
   #       # H[ind_i_p, ind_i_p] = (J_mid'*J_mid + J_next'*J_next).*quad_weight;
   #       H[ind_i_p, ind_i_p] = ((J_mid'*J_mid) + (J_next'*J_next)).*quad_weight;
   #    end
   #    # Off diagonal blocks
   #    for i_p = 1:p
   #       i_p_next = mod1(i_p+1, p);
   #       ind_i_p = 2+(i_p-1)*N:1+i_p*N
   #       ind_i_p_next = 2+(i_p_next-1)*N:1+i_p_next*N
   #
   #       J_mid_next = reshape((i_p_next==1) ? J11 : Jii, 2Nθ, N);
   #       J_next = reshape(J_off_diag[:,:,:,i_p], 2Nθ, N);
   #
   #       H[ind_i_p, ind_i_p_next] = J_mid_next'*J_next .* quad_weight;
   #       H[ind_i_p_next, ind_i_p] = H[ind_i_p, ind_i_p_next]';
   #    end
   # end
end

function full_J(J10, J11, Jii, J_off_diag, z)
   p = get_p(z);
   Nθ = size(J11,2);
   N = get_N(z);
   J = zeros(p*2Nθ, 1+p*N);
   J[1:2Nθ, 1] = J10[:];
   J[1:2Nθ, 2:1+N] = reshape(J11, 2Nθ, N);

   if p > 1
      for i_p = 2:p
         ind_θ = 1+(i_p-1)*2Nθ:i_p*2Nθ
         ind_a_prev = 2 + (i_p-2)*N:1 + (i_p-1)*N
         ind_a = 2 + (i_p-1)*N:1 + i_p*N
         J[ind_θ, ind_a] = reshape(Jii, 2Nθ, N);
         J[ind_θ, ind_a_prev] = reshape(J_off_diag[:,:,:,i_p-1], 2Nθ, N);
      end
      J[1:2Nθ, end-N+1:end] = reshape(J_off_diag[:,:,:,end], 2Nθ, N);
   end

   return J
end


# Gets w = Jᵀv
function J_left_multiply(v::AbstractVector, J10, J11, Jii, J_off_diag, z)
   p = get_p(z);
   N = get_N(z);
   Nθ = size(J11, 2);
   quad_weight = 2π/Nθ;

   w = zeros(1 + p*N);

   if p==1
      w[1] = vec(J10)'*v;
      w[2:N+1] = reshape(J11, 2Nθ, N)'*v;
   else
      w[1] = vec(J10)'*v[1:2Nθ];
      for i_p = 1:p
         i_p_next = mod1(i_p+1, p);
         J_mid = reshape((i_p==1) ? J11 : Jii, 2Nθ, N);
         J_next = reshape(J_off_diag[:,:,:,i_p], 2Nθ, N);

         w_ind = (i_p-1)*N+2:i_p*N+1
         v_ind = (i_p-1)*2Nθ+1:i_p*2Nθ;
         v_ind_next = (i_p_next-1)*2Nθ+1:i_p_next*2Nθ;

         w[w_ind] = J_mid'*v[v_ind] + J_next'*v[v_ind_next]
      end
   end

   return quad_weight .* w;
end

# Gets w = Jv
function J_right_multiply(v::AbstractVector, J10, J11, Jii, J_off_diag, z)
   p = get_p(z);
   N = get_N(z);
   Nθ = size(J11, 2);

   w = zeros(p*2Nθ);

   if p==1
      w[:] = reshape(J10, 2Nθ, 1) * v[[1]] + reshape(J11, 2Nθ, N) * v[2:end];
   else
      w[1:2Nθ] = reshape(J10, 2Nθ, 1) * v[[1]];
      for i_p = 1:p
         ind_w = 1+(i_p-1)*2Nθ:i_p*2Nθ;
         ind_a = 2+(i_p-1)*N:1+i_p*N;
         i_p_next = mod1(i_p+1, p)
         ind_w_next = 1+(i_p_next-1)*2Nθ:i_p_next*2Nθ;

         J_mid = reshape((i_p==1) ? J11 : Jii, 2Nθ, N);
         J_next = reshape(J_off_diag[:,:,:,i_p], 2Nθ, N);
         w[ind_w] += J_mid*v[ind_a]
         w[ind_w_next] += J_next*v[ind_a];
      end
   end

   return w;
end

"""
    get_circle_residual(F::Function, z::InvariantCircle, Nθ::Integer)

Get the KAM-like residual of a chain of p invariant circles.\\
>   Rᵢ₁ = z₁(θᵢ+τ) - F(z_p(θᵢ))\\
>   Rᵢⱼ = zⱼ(θᵢ) - F(zⱼ₋₁(θᵢ)) for 2<=j<=p

Arguments:
- `F`: Symplectic map
- `z`: Invariant circle
- `Nθ`: Number of θ points where the residual is taken
"""
function get_circle_residual(F::Function, z::InvariantCircle, Nθ::Integer)
   # Get parameters
   p  = get_p(z);
   N = get_N(z);
   id = Matrix(I*1.0, 2, 2);

   R = zeros(2, Nθ, p);
   Φ = get_Φ(z, Nθ);
   Φτ = get_Φ(z, Nθ; τ=get_τ(z));

   # Initialize
   x = zeros(2, Nθ, p);
   grid_eval!(x, z, Φ);

   # Fill residual with evaluations of z
   x_τ = zeros(2, Nθ)
   grid_eval!(x_τ, z, Φτ; i_circle=1);
   R[:, :, 1] = x_τ;
   if p>1
      R[:, :, 2:end] = x[:, :, 2:end]
   end


   # variation of F(C(i-1)) and residual
   for i_p = 1:p
      for i_θ = 1:Nθ
         Fi = F(x[:, i_θ, i_p])

         # Adjust residual
         i_R = mod1(i_p+1, p);
         R[:, i_θ, i_R] += -Fi;
      end
   end

   return R;
end

"""
    gn_circle(FJ::Function, z::InvariantCircle; Nθ::Integer=0,
              maxiter::Integer=10, rtol::Number=1e-8, verbose::Bool=false,
              monitor::Function=(z, rnorm) -> nothing,
              constraint::Function=fixed_θ0_constraint, λ::Number=0)

Find an invariant circle using Gauss-Newton with linesearch. Implemented with
dense linear algebra (no FFTs yet).

Arguments:
- `z::InvariantCircle`: An initial connecting orbit guess, see
  `linear_initial_connecting_orbit`
- `FJ::Function`: Function defined on R² with signature
  `F(x), J(x) = FJ(x)` where `F` is the symplectic map and `J = dF/dx`
- `ends::AbstractArray`: A 2p × 2 matrix containing the end periodic orbits
  `[x00, x10; F(x00), F(x10); ... ; F^(p-1)(x00), F^(p-1)(x10)]`
- `Ns = 100`: Number of quadrature nodes for the least squares
- `maxiters = 10`: Maximum number of Gauss Newton steps
- `rtol = 1e-8`: Stopping tolerance
"""
function gn_circle(FJ::Function, z::InvariantCircle; Nθ::Integer=0,
                   maxiter::Integer=10, rtol::Number=1e-8, verbose::Bool=false,
                   monitor::Function=(z, rnorm) -> nothing,
                   constraint::Function=fixed_θ0_constraint,
                   λ::Number=0)
   # Indexing and allocation
   p = get_p(z);
   N = get_N(z);

   R = zeros(2, Nθ, p);
   J10 = zeros(2, Nθ)
   J11 = zeros(2, Nθ, N);
   Jii = zeros(2, Nθ, N);
   J_off_diag = zeros(2, Nθ, N, p);
   Φ = get_Φ(z, Nθ);
   H = zeros(3 + p*N, 3+p*N);

   # Constraints
   c = constraint(z);
   H[1:p*N+1, end-1:end] = c';
   H[end-1:end, 1:p*N+1] = c;
   μ = zeros(2);
   μ0 = c[:, 1]*get_τ(z) + c[:, 2:end]*vec(get_a(z));


   # Iterate
   rnorm = 0.0
   rhs = zeros(3+p*N)

   Φτ = get_Φ(z, Nθ; τ=get_τ(z));
   try
      resid!(R, J10, J11, Jii, J_off_diag, FJ, z, Φ, Φτ)
   catch
      return maxiter+1
   end

   hessian!(H, J10, J11, Jii, J_off_diag, z)
   rhs[1:1+p*N] = - J_left_multiply(vec(R), J10, J11, Jii, J_off_diag, z);
   rhs[end-1:end] .= μ0 - ( c[:, 1]*get_τ(z) + c[:, 2:end]*vec(get_a(z)) );

   rnorm = norm(rhs)/average_radius(z); # TODO: Is this the right relative error?
   normk = 2*rnorm;

   if verbose
      println("Entering gn_circle with rnorm=$rnorm")
   end

   for k = 1:maxiter
      update = (H + λ*Diagonal(H)) \ rhs;
      a_prev = get_a(z);
      τ_prev = get_τ(z);
      μ_prev = μ;


      α = 1.0
      kmax = 3;
      for k = 1:kmax
         set_τ!(z, τ_prev + α*update[1]);
         set_a!(z, a_prev + α .* reshape(update[2:1+p*N], N, p))
         μ = μ_prev + update[end-1:end]

         Φτ = get_Φ(z, Nθ; τ=get_τ(z));
         try
            resid!(R, J10, J11, Jii, J_off_diag, FJ, z, Φ, Φτ)
         catch
            @warn "Resid! failed, probably due to FJ"
            return maxiter+1
         end

         # TODO: Make an actual solver here. This can be really fast,
         # so that I can feel justified in actually writing the complicated
         # version of the algorithm
         hessian!(H, J10, J11, Jii, J_off_diag, z)
         rhs[1:1+p*N] = - J_left_multiply(vec(R), J10, J11, Jii, J_off_diag, z);
         rhs[end-1:end] .= μ0 - ( c[:, 1]*get_τ(z) + c[:, 2:end]*vec(get_a(z)) );
         normk = norm(rhs)/average_radius(z);

         if (normk < rnorm)
            break
         end
         α /= 2
      end
      rnorm = normk

      monitor(z, rnorm)
      if verbose
         println("gn_circle k=$k, rnorm=$rnorm, α=$α, z(0,1) = $(z(0,1))")
         println("   rhs[1] = $(rhs[1]),
   norm(rhs[2:1+p*N]) = $(norm(rhs[2:1+p*N])),
   norm(rhs(end-1:end)) = $(norm(rhs[end-1:end]))")
      end
      if rnorm < rtol
         if verbose
             println("gn_circle converged in $k steps")
         end
         monitor(z, rnorm)
         return k
      end
   end


   monitor(z, rnorm)
   if verbose
      println("Warning: Did not converge in $(maxiter) steps: $(rnorm) > $(rtol)")
   end
   return maxiter+1
end
