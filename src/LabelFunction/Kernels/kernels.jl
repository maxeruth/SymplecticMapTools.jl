function distance_matrix(x::AbstractVector, y::AbstractVector)
    return (x .- y').^2
end

function distance_matrix(x, y)
    d  = size(x, 1);
    Nx = size(x, 2);
    Ny = size(y, 2);

    K = zeros(Nx, Ny);
    for ii = 1:d
        K .+= (x[ii, :] .- y[ii, :]') .^ 2;
    end

    return K
end

include("./squared_exponential.jl")
include("./multiquadric.jl")
include("./fourier_se.jl")
include("./fourier.jl")

"""
    KernelLabel

Constructors:
- `KernelLabel(x::AbstractArray, c::AbstractVector, σ::Number;`
                `kernel::Symbol=:SquaredExponential)`
- `KernelLabel(x::AbstractArray, c::AbstractVector, σ::AbstractArray;`
                `kernel::Symbol=:SquaredExponential)`

An approximately invariant kernel label function.

Constructor elements:
- `x`: The `d` × `2N` array of interpolating points, where
  `x[:, n+N] = F(x[:, n])` for `n<=N` and some symplectic map `F`.
- `c`: The length `2N` vector of coefficients of the kernel function
- `σ`: The length scale of the kernel (if a vector, length scale is different in
  different directions)
- `kernel`: The type of kernel used

The current supported kernels are:
- `:SquaredExponential`: The squared exponential kernel
    K(x, y) = exp(-|x-y|²)
- `:InverseMultiquadric`: The β=1 Inverse Multiquadric kernel
    K(x, y) = (1+|x-y|²)^(-1/2)
- `:FourierSE`: The Fourier × Cartesian squared exponential kernel
    K(x, y) = exp(-sin²(x₁-y₁) - (x₂-y₂)²)
- `:Fourier`: The Fourier × Fourier squared exponential kernel
    K(x, y) = exp(-sin²(x₁-y₁) - sin²(x₂-y₂))
"""
struct KernelLabel
    x::AbstractArray;  # Interpolating points
    c::AbstractVector; # Interpolating weights
    σ::AbstractArray;  # Characteristic length scale
    kernel::Symbol;    # Type of kernel used

    function KernelLabel(x::AbstractArray, c::AbstractVector, σ::AbstractArray;
                            kernel::Symbol=:SquaredExponential)
        @assert kernel ∈ (:SquaredExponential, :InverseMultiquadric, :FourierSE,
                          :Fourier)
        new(x, c, Diagonal(σ), kernel);
    end
end

function KernelLabel(x::AbstractArray, c::AbstractVector, σ::Number;
                        kernel::Symbol=:SquaredExponential)
    #
    KernelLabel(x, c, σ .* ones(size(x,1)); kernel)
end

function get_x(k::KernelLabel); return k.x; end
function get_c(k::KernelLabel); return k.c; end
function get_σ(k::KernelLabel); return k.σ; end
function get_kernel(k::KernelLabel); return k.kernel; end
function get_N(k::KernelLabel); return size(k.x,2)÷2; end


function set_x!(k::KernelLabel, x::AbstractArray); k.x[:] = x[:]; end
function set_c!(k::KernelLabel, c::AbstractVector); k.c[:] = c[:]; end

"""
    get_matrix(k::KernelLabel, x::AbstractArray)

Get the kernel matrix `K[i,j] = K(k.x[:,i], x[:,j])`
"""
function get_matrix(k::KernelLabel, x::AbstractArray)
    σ = get_σ(k);

    if get_kernel(k) == :SquaredExponential
        return SE_matrix(σ\x, σ\get_x(k))

    elseif get_kernel(k) == :InverseMultiquadric
        return inverse_multiquadric_matrix(σ\x, σ\get_x(k))

    elseif get_kernel(k) == :FourierSE
        return fourier_SE_matrix(x, get_x(k), σ)

    elseif get_kernel(k) == :Fourier
        return fourier_matrix(x, get_x(k), σ)

    end
end

"""
    evaluate(k::KernelLabel, x::AbstractArray)

Evaluate the kernel matrix at the columns of x
"""
function evaluate(k::KernelLabel, x::AbstractArray)
    return get_matrix(k, x) * get_c(k); # Default routine
end

"""
    (k::KernelLabel)(x::AbstractArray)

Wrapper for `evaluate(k::KernelLabel, x::AbstractArray)`
"""
function (k::KernelLabel)(x::AbstractArray)
    evaluate(k, x)
end
