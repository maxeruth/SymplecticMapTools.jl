
function SE_matrix(x, y)
    K = distance_matrix(x, y)
    K .= exp.(-K/2.)
    # K .= exp.(-K/2.)

    return K
end

# function SE_matrix(x, y)
#     d  = size(x, 1);
#     Nx = size(x, 2);
#     Ny = size(y, 2);
#     distances = zeros(d, Nx, Ny)
#
#     K = zeros(Nx, Ny);
#     for ii = 1:Nx, jj = 1:Ny
#         distances[:, ii, jj] = x[:, ii] - y[:, jj]
#     end
#     f_K = (v) -> exp(-v'*v/2);
#     K = [f_K(distances[:, ii, jj]) for ii = 1:Nx, jj = 1:Ny]
#     return Symmetric(K)
# end

function SE_matrix_and_derivs(x, y)
    d  = size(x, 1);
    Nx = size(x, 2);
    Ny = size(y, 2);
    distances = zeros(d, Nx, Ny)

    K = zeros(Nx, Ny);
    dK = zeros(d, Nx, Ny)
    for ii = 1:Nx, jj = 1:Ny
        v = x[:, ii] - y[:, jj]
        K[ii, jj] = -v'*v / 2
        dK[:, ii, jj] = -v .* K[ii, jj]
    end

    return K, dK
end
