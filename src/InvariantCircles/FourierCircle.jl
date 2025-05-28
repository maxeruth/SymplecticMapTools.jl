"""
    FourierCircle <: InvariantCircle

Constructors:
- `FourierCircle(Na::Integer; a::AbstractArray=[], τ::Number=0., p::Integer=1)`

Fourier invariant circle data structure (see InvariantCircle), used to compute a
circle z that is invariant of a map Fᵖ(z), where p is the period. This is done
by storing all of the circles Fⁱ(z) for 0<=i<p. The circles are stored as
Fourier series, i.e.
   z(θ) = a₀ + ∑₁ᴺᵃ Aⱼ [cos(jθ), sin(jθ)]ᵀ
The coefficients of the invariant circle can be accessed and set in chunks via
array indexing. That is,
  You can get and set a₀ of the ith circle using `z[0, i]`
  You can get and set Aⱼ of the ith circle using `z[j, i]`
"""
struct FourierCircle <: InvariantCircle
   a::AbstractArray; # Representation of all of the circles
   τ::AbstractArray; # Rotation number of the circle
   sz::Tuple;        # The size (Na, p)

   function FourierCircle(Na::Integer; a::AbstractArray=[], τ::Number=0., p::Integer=1)
      if (a == [])
         a = zeros((2 + 4*Na), p);
      end

      return new(a, [1.0*τ], (Na, p));
   end
end

function FourierCircle(tor::FourierTorus)
   sz = size(tor.a)
   @assert length(sz) == 3 # Needs to be a 1D torus
   @assert sz[1] == 2      # Needs to be a 1D torus in a 2D domain
   p = sz[2]
   Na = sz[3]÷2
   
   a = zeros(2,2Na+1,p)
   for ii = 1:p
      for jj = 1:2
         a[jj,1,ii] = real.(tor.a[jj,ii,Na+1])
         a[jj,2:2:end,ii] = 2 .* real.(tor.a[jj,ii,Na+2:end])
         a[jj,3:2:end,ii] = -2 .* imag.(tor.a[jj,ii,Na+2:end])
      end
   end

   FourierCircle(Na; a=reshape(a,4Na+2,p), p=p, τ=tor.τ[1])
end

## Getters and setters
function Base.size(z::FourierCircle)
   return z.sz;
end

"""
    get_Na(z::FourierCircle)

Get number of Fourier modes
"""
function get_Na(z::FourierCircle)
   return z.sz[1];
end

"""
    get_N(z::FourierCircle)

Get number of unknown parameters per circle in island, excluding τ
"""
function get_N(z::FourierCircle)
   Na = get_Na(z);
   return 4Na + 2;
end

"""
   `get_p(z::FourierCircle)`

Get number of islands
"""
function get_p(z::FourierCircle)
   return z.sz[2]
end

"""
   `Base.length(z::FourierCircle)`

Get total number of parameters, excluding τ
"""
function Base.length(z::FourierCircle)
   return length(z.a);
end

"""
    get_a0(z::FourierCircle; i_circle::Integer=1)

Get constant term of Fourier series of the `i_circle`th circle in island chain
"""
function get_a0(z::FourierCircle; i_circle::Integer=1)
   return z.a[1:2, i_circle];
end

"""
    set_a0!(z::FourierCircle, a0::AbstractArray; i_circle::Integer=1)

Set constant term of Fourier series of the `i_circle`th circle in island chain
"""
function set_a0!(z::FourierCircle, a0::AbstractArray; i_circle::Integer=1)
   z.a[1:2, i_circle] = a0[:];
end

"""
    get_Am(z::FourierCircle, m::Integer; i_circle::Integer=1)

Get `m`th Fourier coefficients of the `i_circle`th circle in island chain
"""
function get_Am(z::FourierCircle, m::Integer; i_circle::Integer=1)
   return z.a[4*m-2 .+ (1:4), i_circle];
end

"""
    set_Am!(z::FourierCircle, m::Integer, a::AbstractArray; i_circle::Integer=1)

Set `m`th Fourier coefficients of the `i_circle`th circle in island chain
"""
function set_Am!(z::FourierCircle, m::Integer, a::AbstractArray; i_circle::Integer=1)
   z.a[4*m-2 .+ (1:4), i_circle] = a[:];
end

"""
    Base.getindex(z::FourierCircle, i_A::Integer, i_circle::Integer)

Get the coefficients of z. `z[0, i_circle]` gets the constant term.
`z[i_A, i_circle]` gets higher coefficients
"""
function Base.getindex(z::FourierCircle, i_A::Integer, i_circle::Integer)
   if i_A == 0
      return get_a0(z; i_circle)
   else
      return get_Am(z, i_A; i_circle)
   end
end

"""
    Base.setindex!(z::FourierCircle, X::AbstractArray, i_A::Integer,
                   i_circle::Integer)

Set the coefficients of z. `z[0, i_circle]` sets the constant term.
`z[i_A, i_circle]` sets higher coefficients
"""
function Base.setindex!(z::FourierCircle, X::AbstractArray, i_A::Integer, i_circle::Integer)
   if i_A == 0
      set_a0!(z, X; i_circle);
   else
      set_Am!(z, i_A, X; i_circle)
   end
end

function get_a(z::FourierCircle; i_circle::Integer=0)
   if i_circle == 0
      return deepcopy(z.a);
   else
      return z.a[:, i_circle]
   end
end

function set_a!(z::FourierCircle, a::AbstractArray; i_circle::Integer=0)
   if i_circle == 0
      z.a[:] = a[:];
   else
      z.a[:, i_circle] = a[:];
   end
end

"""
    get_τ(z::FourierCircle)

Get the rotation number (the circle is 2π periodic)
"""
function get_τ(z::FourierCircle)
   return z.τ[1];
end

"""
    set_τ!(z::FourierCircle, τ::Number)

Set the rotation number (the circle is 2π periodic)
"""
function set_τ!(z::FourierCircle, τ::Number)
   z.τ[1] = τ;
end

"""
    average_radius(z::FourierCircle)

Norm defined as the L2 norm of the radius, i.e.
‖z‖² = 1/2πp ∫ ||z(θ) - a₀||^2 dθ
     = 1/2p ∑ₖ Tr(AₖᵀAₖ)
If `i_circle = 0`, takes the norm over all circles. If `i_circle != 0`, finds
the average radius of a single circle.
"""
function average_radius(z::FourierCircle; i_circle::Integer = 0)
   ps = (i_circle == 0) ? (1:get_p(z)) : [i_circle]
   norm_sq = 0.0;
   for i_A = 1:get_Na(z), i_p = ps
      A_i = reshape(z[i_A, i_p], 2, 2);
      norm_sq += tr(A_i'*A_i)
   end
   return sqrt(norm_sq/(2*sum(ps)));
end

"""
    evaluate(z::FourierCircle, θ::AbstractVector{T}; i_circle::Integer=1) where {T}

Evaluate the `i_circle`th circle in the chain at a vector of points θ
"""
function evaluate(z::FourierCircle, θ::AbstractVector{T}; i_circle::Integer=1) where {T}
   Nθ = length(θ);
   x = zeros(T, 2, Nθ);
   ϕ = zeros(T, 2, Nθ);

   for jj = 1:Nθ
      x[:,jj] = z[0, i_circle];
   end

   for i_A = 1:get_Na(z)
      ϕ[1,:] = cos.(i_A .* θ);
      ϕ[2,:] = sin.(i_A .* θ);
      x += reshape(z[i_A, i_circle], 2, 2) * ϕ;
   end

   return x
end

"""
    deval(z::FourierCircle, θ::AbstractVector; i_circle::Integer=1)

Evaluate the derivative of the `i_circle`th circle in the chain at a vector of
points `θ`
"""
function deval(z::FourierCircle, θ::AbstractVector; i_circle::Integer=1)
   Nθ = length(θ);
   dx = zeros(2, Nθ);
   dϕ = zeros(2, Nθ);

   for i_A = 1:get_Na(z)
      dϕ[1,:] = -sin.(i_A.*θ);
      dϕ[2,:] =  cos.(i_A.*θ);
      dx += reshape(i_A.*z[i_A, i_circle], 2, 2) * dϕ;
   end

   return dx
end

function Base.similar(z::FourierCircle)
   FourierCircle(get_Na(z); p=get_p(z))
end

"""
    deriv(z::FourierCircle)

Return the derivative of the FourierCircle, wrapped in a circle object
"""
function deriv(z::FourierCircle)
   R = [0. -1.; 1. 0.]

   zp = similar(z)
   set_τ!(zp, get_τ(z))

   for i_A = 1:get_Na(z)
      for i_circle = 1:get_p(z)
         zi = reshape(z[i_A, i_circle], 2, 2)
         zp[i_A, i_circle] = reshape(i_A.*zi*R, 4)
      end
   end

   zp
end

## Accelarated routines
# Routines that work on a uniform grid, using a matrix Φ which contains all
# of the Fourier terms
function get_Φ(z::FourierCircle, Nθ::Integer; τ::Number=0.0)
   θ = ((0:Nθ-1) .* (2π/Nθ)) .+ τ;
   Φ = zeros(Nθ, get_Na(z)*2 + 1);
   Φ[:, 1] .= 1.;
   for iA = 1:get_Na(z)
      Φ[:, 2*iA]   = cos.(iA.*θ);
      Φ[:, 2*iA+1] = sin.(iA.*θ);
   end

   return Φ;
end

# x should be size 2×Nθ×p
function grid_eval!(x::AbstractArray, z::FourierCircle, Φ::AbstractArray; i_circle=0)
   Nθ = size(Φ, 1);

   if i_circle==0
      for i_p = 1:get_p(z)
         x[:, :, i_p] = z[0, i_p]*Φ[:, [1]]';
      end
   else
      x[:,:] = z[0, i_circle]*Φ[:, [1]]';
   end

   for i_A = 1:get_Na(z)
      if i_circle==0
         for i_p = 1:get_p(z)
            x[:, :, i_p] += reshape(z[i_A, i_p], 2, 2) * Φ[:, 2i_A:2i_A+1]';
         end
      else
         x[:, :] += reshape(z[i_A, i_circle], 2, 2) * Φ[:, 2i_A:2i_A+1]';
      end
   end
end

# dx should be size 2×Nθ×p
function grid_deval!(dx::AbstractArray, z::FourierCircle, Φ::AbstractArray; i_circle=0)
   Nθ = size(Φ, 2);
   d_mat = [0. -1.; 1. 0.]
   dx[:] .= 0.;
   for i_A = 1:get_Na(z)
      if i_circle == 0
         for i_p = 1:get_p(z)
            dx[:, :, i_p] += (i_A.*reshape(z[i_A, i_p], 2, 2)*d_mat) * Φ[:, 2*i_A:2*i_A+1]';
         end
      else
         dx[:, :] += (i_A.*reshape(z[i_A, i_circle], 2, 2)*d_mat) * Φ[:, 2*i_A:2*i_A+1]';
      end
   end
end


"""
    area(z::FourierCircle; i_circle=1)

Get the area of the circle.
"""
function area(z::FourierCircle; i_circle=1)
   Na = get_Na(z);
   A = 0;
   for i_A = 1:Na
      Am = reshape(get_Am(z, i_A; i_circle), 2, 2)
      A += i_A * π * det(Am);
   end
   return A
end


## Constraints
function fixed_θ0_constraint(z::FourierCircle)
   Na = get_Na(z);
   p = get_p(z);
   c = zeros(2, 1 + p*(4*Na+2))

   id = Matrix(1.0*I, 2, 2);

   # Fix z(0)
   c[:, 2:3] = id;
   for ii = 1:Na
     c[:,4*ii:4*ii+1] = id;
   end

   return c
end

# Average radius is ‖z‖² = ∑‖Aⱼ‖²,
# so d‖z‖/dA = A/‖z‖
function average_radius_constraint(z::FourierCircle)
   Na = get_Na(z);
   N = get_N(z);
   p = get_p(z);
   c = zeros(1, 1+p*N);

   znorm = norm(z);
   for i_A = 1:Na
      c[1, 4i_A : 4i_A+3] = vec(z[i_A, 1])./znorm;
   end

   return c;
end

function area_constraint(z::FourierCircle)
   Na = get_Na(z);
   N = get_N(z);
   p = get_p(z);
   c = zeros(1, 1+p*N);

   D = Diagonal([1,-1,-1,1]); # A = ∑ᵢ Aᵢᵀ D Aᵢ[4:-1:1]
   for i_A = 1:Na
      zi = z[i_A, 1]
      c[1, 4i_A : 4i_A+3] = i_A * π * D * zi[4:-1:1]
   end

   return c
end


"""
    circle_linear!(z::FourierCircle, FJ::Function, x::AbstractArray, h::Number)

Fill `z`'s Fourier coefficients using a linearization about a stable periodic
orbit. Useful for seeding continuation.

Arguments:
- `z`: An invariant circle
- `FJ`: A function that returns the map and its derivative `F, dFdx = FJ(x)`
- `x`: A periodic orbit of `F`
- `h`: The average radius of the first linearized invariant circle
"""
function circle_linear!(z::FourierCircle, FJ::Function, x::AbstractArray,
                        h::Number)
    p = get_p(z);
    if (size(x, 2) != p)
        println("circle_linear!: z should be the same period as x!")
        return
    end

    for i_p = 1:p
        z[0, i_p] = x[:, i_p];
    end

    # Get Jacobians
    Js = zeros(2, 2, p)
    J_full = Matrix(1.0*I, 2, 2);
    for i_p = 1:p
        Js[:,:,i_p] = FJ(x[:, i_p])[2];
        Js[:,:,i_p] = Js[:,:,i_p]
        J_full = Js[:,:,i_p] * J_full;
    end


    eigenJ = eigen(J_full)

    # Choose eigenvector for positive orientation
    if det([real(eigenJ.vectors[:,1]) imag(eigenJ.vectors[:,1])]) > 0
        v = eigenJ.vectors[:,1]
        λ = eigenJ.values[1]
    else
        v = eigenJ.vectors[:,2]
        λ = eigenJ.values[2]
    end

    # Rotate and scale so z(0) = e1
    ur = real(v)
    ui = imag(v)
    r = sqrt(ur[2]^2 + ui[2]^2)
    c = ui[2]/r
    s = -ur[2]/r
    A1 = [ur ui] * [c -s ; s c]
    A1[2,1] = 0.0

    z[1,1] = A1;
    z[1,1] = A1 .* h/average_radius(z; i_circle=1)

    # Return rotation angle and normalized A1 coefficient matrix
    τ = atan(-imag(λ), real(λ))
    set_τ!(z, τ);

    # Fill the rest of the invariant circles
    for i_p = 2:p
        z[1, i_p] = Js[:, :, i_p-1]*reshape(z[1, i_p-1], 2, 2)#/(det(J_full)^(1/(2p)));
    end

    Na = get_Na(z);
    for i_A=2:Na, i_p = 1:p
        z[i_A, i_p] = zeros(4);
    end

    return τ, A1
end
