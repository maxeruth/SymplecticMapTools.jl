function multiquadric_matrix(x, y; β = 1)
    K = distance_matrix(x, y);
    if β == 1
        K .= sqrt.(1. .+ K)
    else
        K .= (1. .+ K) .^ (β/2);
    end

    return K
end

function inverse_multiquadric_matrix(x, y; β = 1)
    return 1 ./ multiquadric_matrix(x, y; β)
end
