function c_to_y(c, K, dK)
    y = K*c;
    dy = dK*c;
    return y, dy;
end

function reguralization_info!(resid, A11, c, y, η, N; get_g=true, get_h=true)
    ind_y = 1:3:6N;

    get_g && (resid[ind_y] = resid[ind_y] + η.*c);
    get_h && (A11[ind_y, :] = A11[ind_y, :] + η*I);

    return (η/2)*y'*c;
end

function constraint_info!(resid, G, y, λ, N; get_g=true, get_h=true)
    Ny = 2N;
    Ng = 6N;

    ind_y = 1:3:Ng;
    ind_λ = Ng+1:Ng+N;

    if get_g
        resid[ind_y] = resid[ind_y] + G*λ;
        resid[ind_λ] = resid[ind_λ] + G'*y;
    end

    return  λ'*G'*y
end

# If s = 0, do nothing. If s = 1, the full problem is L1_info!
function L1_info!(resid, A11, K, y, ytil, s, N; get_g=true, get_h=true)
    ind_y = 1:3:6N;

    get_g && (resid[ind_y] = resid[ind_y] + s.*(y-ytil));
    get_h && (A11[ind_y,:] = A11[ind_y, :] + s.*K);

    return (s/2) * (y-ytil)'*(y-ytil);
end

function L2_info!(resid, A11, K, dK, ∇y, s, N; get_g=true, get_h=true)
    Ny = 2N;
    Ng = 6N;

    ind_y=falses(Ng+N);
    ind_y[1:3:Ng] .= true;
    ind_∇y = .!ind_y;
    ind_∇y[Ng+1:end] .= false;

    ∇y = reshape(∇y, 2, Ny);
    ∇y2m1 = [∇y[:, ii]'*∇y[:,ii] - 1 for ii = 1:Ny]; # "|∇yᵢ|² - 1"

    get_g && (resid[ind_∇y] = resid[ind_∇y] + s .* vec(∇y * Diagonal(∇y2m1)));

    if get_h
        dK = reshape(dK, 2, Ny, Ny);
        for ii = 1:Ny
            ind_ii = (ii-1)*3 .+ (2:3);
            K_ii = 2*(∇y[:,ii]*∇y[:,ii]') + ∇y2m1[ii]*I;
            A11[ind_ii, :] = A11[ind_ii, :] + (s.*K_ii) * dK[:, ii, :]
        end
    end

    return (s/4) * ∇y2m1'*∇y2m1;
end

function continuation_info!(resid, A11, G, K, dK, c, λ, y, ∇y, ytil, s, η, N; get_g=true, get_h=true);
    L = 0.;
    L += reguralization_info!(resid, A11, c, y, η, N; get_g, get_h)
    L += constraint_info!(resid, G, y, λ, N; get_g, get_h)
    L += L1_info!(resid, A11, K, y, ytil, 1-s, N; get_g, get_h)
    L += L2_info!(resid, A11, K, dK, ∇y, s, N; get_g, get_h)
    return L
end

function kernel_solve(resid, A11_in, A12_in, A21, QK)
    Ng, Ny = size(A11_in);
    Nλ = size(A12_in, 2);

    # println("size(A11) = $(size(A11))")
    # println("size(QK) = $(size(QK))")
    A11 = QK'*A11_in
    E1 = QK'*resid[1:Ng];
    E2 = resid[Ng+1:end]
    A12 = QK[1:3:end, :]'*A12_in;
    # println("size(A11) = $(size(A11))")

    QR = qr(A11);
    Q = QR.Q;
    R = QR.R;
    W = Q * (R' \ A21');


    rhs = - W'*E1 + E2;
    lhs = W'*A12;

    Δλ = lhs \ rhs;

    rhs_d = -E1;
    rhs_d = rhs_d - A12*Δλ;
    Δc = QR \ rhs_d
    # println("cond(R) = $(cond(QR.R))")

    return Δc, Δλ
end

function stacked_K(K, dK)
    Ny = size(K,1);
    fullK = zeros(3Ny, Ny);
    fullK[1:3:end, :] = K;
    fullK[2:3:end, :] = dK[1:2:end];
    fullK[3:3:end, :] = dK[2:2:end];

    return fullK
end

function kernel_Q(K, dK)
    Ny = size(K,1);
    fullK = zeros(3Ny, Ny);
    fullK[1:3:end, :] = K;
    fullK[2:3:end, :] = dK[1:2:end];
    fullK[3:3:end, :] = dK[2:2:end];

    QR = qr(fullK);

    return Matrix(QR.Q)
end

function newton_kernel!(k::KernelLabel, λ, s, η; rtol = 1e-6, N_iter=10, verbose=true)
    Nλ = length(λ);
    Ny = 2Nλ;
    Ng = 3Ny;

    K, dK = matrix_and_derivs(k)
    QK = kernel_Q(K, dK);
    # QK = zeros(Ng, Ny);
    # for ii = 1:Ny
    #     QK[3*(ii-1)+1, ii] = 1.0;
    # end
    x = get_x(k);
    G = constraint_matrix(Nλ)
    ytil = [norm(x[:, ii]) for ii = 1:Ny];

    A11 = zeros(Ng, Ny)
    A12 = G;
    A21 = G'*K
    resid = zeros(Ng + Nλ);
    y, ∇y = c_to_y(c, K, dK)

    get_err = (resid) -> begin
        err_vec = zeros(Ny+Nλ);
        err_vec[1:Ny] = QK'*resid[1:Ng];
        err_vec[Ny+1:end] = resid[Ng+1:end];
        return norm(err_vec)
    end

    L = continuation_info!(resid, A11, G, K, dK, c, λ, y, ∇y, ytil, s, η, Nλ);
    err = get_err(resid);
    println("Entering Newton loop err=$err")

    for ii = 1:N_iter
        Δc, Δλ = kernel_solve(resid, A11, A12, A21, QK)
        c[:] = c + Δc;
        λ[:] = λ + Δλ;
        y, ∇y = c_to_y(c, K, dK)
        resid[:] .= 0.;
        A11[:] .= 0.;

        L = continuation_info!(resid, A11, G, K, dK, c, λ, y, ∇y, ytil, s, η, Nλ);
        err = get_err(resid);

        if verbose
            resid_multiplied = QK'*resid[1:3Ny]
            println("ii=$ii, err=$err, projectec_err = $(norm(resid_multiplied)), λ_err = $(norm(resid[Ng+1:end]))")
        end

        if err < rtol
            set_c!(k, c)
            return L
        end
    end

    println("newton_kernel did not converge after $N_iter iterations")
    set_c!(k, c)
    return -1
end


function opt_kernel!(k::KernelLabel, λ, s, η; verbose=true)
    Nλ = length(λ);
    Ny = 2Nλ;
    Ng = 3Ny;

    K, dK = matrix_and_derivs(k)
    QK = kernel_Q(K, dK);
    # QK = stacked_K (K, dK)
    x = get_x(k);
    G = constraint_matrix(Nλ)
    ytil = [norm(x[:, ii]) for ii = 1:Ny];

    A12 = QK[1:3:end, :]'*G;
    A21 = G'*K

    function f_opt(cλ)
        c = cλ[1:Ny];
        λ = cλ[Ny+1:end]
        y, ∇y = c_to_y(c, K, dK)
        return continuation_info!(nothing, nothing, G, K, dK, c, λ, y, ∇y, ytil, s, η, Nλ, get_g=false, get_h=false);
    end

    function g_opt!(resid, cλ)
        c = cλ[1:Ny];
        λ = cλ[Ny+1:end]
        y, ∇y = c_to_y(c, K, dK)
        resid_long = zeros(Ng + Nλ)
        continuation_info!(resid_long, nothing, G, K, dK, c, λ, y, ∇y, ytil, s, η, Nλ, get_g=true, get_h=false);
        resid[1:Ny] = QK'*resid_long[1:Ng]
        resid[Ny+1:end] = resid_long[Ng+1:end];
    end

    function h_opt!(hess, cλ)
        c = cλ[1:Ny];
        λ = cλ[Ny+1:end]
        y, ∇y = c_to_y(c, K, dK)
        A11 = zeros(Ng, Ny)
        continuation_info!(nothing, A11, G, K, dK, c, λ, y, ∇y, ytil, s, η, Nλ, get_g=false, get_h=true);
        hess[1:Ny, 1:Ny] = QK'*A11;
        hess[Ny+1:end, 1:Ny] = A21;
        hess[1:Ny, Ny+1:end] = A12;
        hess[Ny+1:end, Ny+1:end] .= 0.0;
    end

    cλ0 = zeros(Ny + Nλ);
    cλ0[1:Ny] = get_c(k);
    cλ0[Ny+1:end] = λ;

    # resid1 = zeros(Ny + Nλ);
    # resid2 = copy(resid1);
    # hess1 = zeros(Ny+Nλ, Ny+Nλ);
    #
    # ϵ = 1e-4
    # dcλ = ϵ*randn(Ny+Nλ);
    # cλ2 = cλ0 + dcλ;
    #
    # f1 = f_opt(cλ0);
    # g_opt!(resid1, cλ0);
    # h_opt!(hess1, cλ0);
    #
    # f2 = f_opt(cλ2);
    # g_opt!(resid2, cλ2)
    # println("f2-f1       = $(f2-f1)")
    # println("resid1'*dcλ = $(resid1'*dcλ)")
    # println("resid2-resid1:"); flush(stdout); display(resid2-resid1);
    # println("hess1*dcλ:"); flush(stdout); display(hess1*dcλ);


    res = optimize(f_opt, g_opt!, h_opt!, cλ0)

    cλend = Optim.minimizer(res);
    set_c!(k, cλend[1:Ny]);
    λ[:] = cλend[Ny+1:end];

    return res
end
