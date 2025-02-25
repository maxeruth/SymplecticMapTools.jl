abstract type InvariantTorus end
    
"""
    FourierTorus <: InvariantTorus

Constructors:
- `FourierTorus(d::Integer, Na::Vector{<:Integer}; a::AbstractArray=[], τ::AbstractVector=[0.], p::Integer=1)`

Fourier invariant torus data structure (see InvariantTorus), used to compute a
torus z that is invariant of a map Fᵖ(z), where p is the period. This is done
by storing all of the tori Fⁱ(z) for 0<=i<p. The tori are stored as
Fourier series, i.e.
   z(θ) = a₀ + ∑₁ᴺᵃ Aⱼ [cos(jθ), sin(jθ)]ᵀ
The coefficients of the invariant torus can be accessed and set in chunks via
array indexing. That is,
  You can get and set a₀ of the ith torus using `z[0, i]`
  You can get and set Aⱼ of the ith torus using `z[j, i]`
"""
struct FourierTorus <: InvariantTorus
    a::AbstractArray; # Representation of all of the circles
    τ::AbstractVector; # Rotation number of the circle
    sz::Tuple;        # The size (Na, p)
    
    function FourierTorus(d::Integer, Na::Vector{<:Integer}; a::AbstractArray=[], τ::AbstractVector=[0.], p::Integer=1)
        if (a == [])
            a = zeros(Complex{Float64}, d, p, (1 .+ 2 .* Na)...);
        end

        if τ == [0.]
            τ = zeros(length(Na))
        end
        
        return new(a, τ, (d,p,Na));
    end
end

function evaluate(tor::FourierTorus, theta::AbstractVector)
    d,p,Na = tor.sz
    D = length(Na)

    a = tor.a
    for ii = D:-1:1
        a = reshape(a, length(a)÷(2Na[ii]+1), (2Na[ii]+1))
        F = [exp(im*n*theta[ii]) for n = -Na[ii]:Na[ii]]
        a = a*F
    end
    real.(reshape(a,d,p))
end

"""
    evaluate_on_grid(tor::FourierTorus, thetavecs::AbstractVector)


"""
function evaluate_on_grid(tor::FourierTorus, thetavecs::AbstractVector)
    # Get size
    d,p,Na = tor.sz
    D = length(Na)

    # "Unroll" the Fourier evaluation on each dimension
    a = tor.a
    Nthetas = vcat(length.(thetavecs), 1)
    for ii = D:-1:1
        Ntheta_ii = prod(Nthetas[ii+1:end])
        Na_prev = (ii==1) ? d*p : d*p*prod(2 .* Na[1:ii-1] .+ 1) 
        a = reshape(a, Na_prev, 2Na[ii]+1, Ntheta_ii)
        
        F = [exp(im*n*theta) for n = -Na[ii]:Na[ii], theta in thetavecs[ii]]

        new_a = zeros(Complex{Float64},Na_prev, Nthetas[ii], Ntheta_ii)
        for jj = 1:Ntheta_ii
            new_a[:,:,jj] = a[:,:,jj]*F
        end
        a = new_a
    end
    real.(reshape(a,d,p,Nthetas[1:end-1]...))
end