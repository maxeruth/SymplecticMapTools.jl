function fourier_SE_matrix(x, y, σ)
    σ = diag(σ)

    K = π .* sqrt.(distance_matrix(x[1,:], y[1,:]))
    K = (sin.(K)./(σ[1]*π)).^2;
    K = K + distance_matrix(x[2,:], y[2,:])./(σ[2]^2)

    K .= -K/2.;
    K .= exp.(K)

    return K
end
