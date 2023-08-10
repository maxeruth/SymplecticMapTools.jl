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
derivative (i.e. `F, dF/dx = FJ(x)`) with parameter `k`.
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
    polar_map(;z0 = -0.5)

Returns the polar map h:(θ,z)↦((z-z0)cos(2πθ), (z-z0)sin(2πθ)), as well as
its derivative `HJ`, its inverse `hinv`, and its inverse derivative `HJinv`.
Useful for applying extrapolation methods to maps on 𝕋×ℝ. Default value of
`z0` is useful for the standard map on 𝕋×[0,1] with k=0.7.
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
