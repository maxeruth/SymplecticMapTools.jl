## Standard Map
"""
    standard_map_F(k)

Returns the Chirikov standard map with parameter `k`.
"""
function standard_map_F(k::Number)
    function F(xin)
        x = xin[1]
        y = xin[2]

        yp = y - (k/(2π))*sin(2π*x)
        return [mod(x + yp, 1), yp];
    end

    return F
end

"""
    standard_map_FJ(k)

Returns a function `FJ` that returns the Chirikov standard map and its
derivative (i.e. `F, dFdx = FJ(x)`) with parameter `k`.
"""
function standard_map_FJ(k)
    function FJ(xin)
        x, y = xin;

        g = y - (k/(2π))*sin(2π*x);
        gx = -k*cos(2π*x);
        gy = 1;
        f = mod(x + g, 1);
        fx = 1 + gx;
        fy = gy;

        return [f, g], [fx fy; gx gy]
    end

    return FJ
end


"""
    standard_map_F(K::AbstractMatrix, delta::AbstractVector)

Returns the a higher dimensional standard map with matrix `K` and offset `delta`,i.e.
  y+ = y - (K/(2π))*sin(2π(x+δ))
  x+ = x + y+
where (x,y) take up the (1:d,d+1:2d) entries of the input.
"""
function standard_map_F(K::AbstractMatrix, delta::AbstractVector)
    d = length(delta)
    @assert size(K) == (d,d)
    K2p = K ./ (2π)

    function F(xin)
        x = xin[1:d]
        y = xin[d+1:2d]

        yp = y - K2p*sin.(2π*(x+delta))
        return vcat(mod.(x + yp, 1), yp);
    end

    return F
end

"""
    standard_map_FJ(K::AbstractMatrix, delta::AbstractVector)

Returns a function `FJ` that returns the higher dimensional standard map and its
derivative (i.e. `F, dFdx = FJ(x)`) with parameter `k` 
(see [`standard_map_F(::AbstractMatrix,::AbstractVector)`](@ref)).
"""
function standard_map_FJ(K::AbstractMatrix, delta::AbstractVector)
    d = length(delta)
    @assert size(K) == (d,d)
    K2p = K ./ (2π)
    indx = 1:d
    indy = 1+d:2d

    function FJ(xin)
        x = xin[1:d]
        y = xin[d+1:2d]

        f = zeros(2d)
        df = zeros(2d,2d)
        f[indy] = y - (K * sin.(2π .* (x+delta))) ./ (2π)
        f[indx] = mod.(x + f[indy], 1);
        
        gx = -K*Diagonal(cos.(2π*(x+delta)));
        
        df[indx,indx] = I + gx;
        df[indx,indy] = df[indx,indy] + I
        df[indy,indx] = gx
        df[indy,indy] = df[indy,indy] + I

        return f, df
    end

    return FJ
end


"""
    polar_map(d, z0)

Returns the d-dimensional polar map \\
> `h:(θ,z)->((z-z0)cos(2πθ), (z-z0)sin(2πθ))`\\
as well as its derivative `HJ`, its inverse `hinv`, and its inverse derivative
`HJinv`. Useful for applying extrapolation methods to maps on Tᵈ×Rᵈ. 
"""
function polar_map(d::Integer, z0::AbstractVector)
    ix = 1:d
    iy = d+1:2d
    function h(xin)
        xout = zeros(2d)

        theta = 2π .* xin[ix]
        R = xin[iy] - z0

        xout[ix] = R .* cos.(theta)
        xout[iy] = R .* sin.(theta)
        
        xout
    end

    function HJ(xin)
        xout = zeros(2d)
        dxout = zeros(2d,2d)

        theta = 2π .* xin[ix]
        R = xin[iy] - z0
        
        xout[ix] = R .* cos.(theta)
        xout[iy] = R .* sin.(theta)
        
        dxout[ix,ix] = Diagonal(- 2π .* R .* sin.(theta))
        dxout[iy,ix] = Diagonal(  2π .* R .* cos.(theta))
        dxout[ix,iy] = Diagonal(cos.(theta))
        dxout[iy,iy] = Diagonal(sin.(theta))

        xout, dxout
    end
    
    function hinv(xin)
        xout = zeros(2d)
        
        x = xin[ix]
        y = xin[iy]

        xout[ix] = [mod(atan(y[ii], x[ii])/(2π), 1) for ii = 1:d]
        xout[iy] = sqrt.(x.^2 + y.^2) + z0
        
        xout
    end
    
    function HJinv(x)
        f = hinv(x);
        _, b = HJ(f);
        f, inv(b)
    end

    return h, HJ, hinv, HJinv
end

"""
    polar_map([d];z0 = -0.5)

Returns the polar map \\
> `h:(θ,z)->((z-z0)cos(2πθ), (z-z0)sin(2πθ))`\\
as well as its derivative `HJ`, its inverse `hinv`, and its inverse derivative
`HJinv`. Useful for applying extrapolation methods to maps on T×R. Default
value of `z0` is useful for the standard map on T×[0,1] with k=0.7.
"""
function polar_map(;z0 = -0.5)
    h = (x) -> [(x[2]-z0)*cos(2π*x[1]), (x[2]-z0)*sin(2π*x[1])]
    HJ = (x) -> (h(x), [-2π*(x[2]-z0)*sin(2π*x[1]) cos(2π*x[1]);
                         2π*(x[2]-z0)*cos(2π*x[1]) sin(2π*x[1])]);
    hinv = (x) -> [mod(atan(x[2], x[1])/(2π), 1), sqrt(x[1]^2 + x[2]^2) + z0];
    HJinv = (x) -> begin
        f = hinv(x);
        _, b = HJ(f);
        f, inv(b)
    end

    return h, HJ, hinv, HJinv
end