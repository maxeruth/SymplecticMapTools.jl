function fourier_SE_matrix(x, y, σ)
    d = size(x, 1);
    Nx = size(x, 2);
    Ny = size(y, 2);
    @assert size(y, 1) == d;
    @assert d == 2

    σ = diag(σ)
    @assert length(σ) == 2
    σ1 = σ[1];
    σ2 = σ[2];

    K = zeros(Nx, Ny)

    @inbounds for jj = 1:Ny
        y1 = y[1,jj];
        y2 = y[2,jj];
        @inbounds for ii = 1:Nx
            x1 = x[1,ii];
            x2 = x[2,ii];

            d1 = (sin(π*abs(x1-y1))/(π*σ1))^2
            d2 = ((x2-y2)/σ2)^2
            K[ii, jj] = exp(-(d1+d2)/2)
        end
    end

    return K
end
