function fourier_matrix(x, y, σ)
    σ = diag(σ)

    d1 = π .* sqrt.(distance_matrix(x[1,:], y[1,:]))
    d2 = π .* sqrt.(distance_matrix(x[2,:], y[2,:]))
    K = (sin.(d1)./(σ[1]*π)).^2 + (sin.(d2)./(σ[2]*π)).^2;

    K .= -K/2.;
    K .= exp.(K)

    return K
end
