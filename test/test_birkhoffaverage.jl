@testset "BirkhoffAveraging.jl" begin
    k_sm = 0.2;
    x0 = [0.49, 0.5]

    F = standard_map_F(k_sm)
    h, HJ, hinv, HJinv = polar_map()

    # vector_mpe_backslash, vector_mpe_iterative, ContFrac, big_cont_frac_eval,
    #        partial_frac, denoms, big_partial_frac, big_denoms, wba_weight,
    #        birkhoff_extrapolation, adaptive_birkhoff_extrapolation, sum_stats,
    #        get_sum_ave, get_circle_info

    Kinit = 10
    Kstride = 50
    Kmax = 110

    ## Island Test
    sol = adaptive_birkhoff_extrapolation(h, F, x0; Kinit, Kmax, Kstride, Nfactor=1)
    
    # Birkhoff average should match the rotation number
    @test abs(dot(sol.xs[2,1:2sol.K+1],sol.c) - 0.5) < 1e-8

    # Test the island number is correct
    get_w0!(sol, 1)
    @test denominator(sol.w0[1]==2) 

    # Test the torus is correct
    adaptive_get_torus!(sol)
    Ntheta = 10
    @test norm(get_circle_residual((y) -> h(F(hinv(y))), sol.tor, Ntheta)) < 1e-10

    ## Dimension tests
    w0s = mod.([sqrt(5)/2-1, 2-sqrt(3), sqrt(2)-1], 1.)
    k_sms = [zeros(1,1), zeros(2,2), zeros(3,3)]
    deltas = [zeros(1), zeros(2), zeros(3)]
    for dim = 1:3
        F = standard_map_F(k_sms[dim],deltas[dim])
        z0 = ones(dim) .* -0.5
        h, HJ, hinv, HJinv = polar_map(dim, z0)

        x0 = zeros(2dim)
        x0[dim+1:end] = w0s[1:dim]
        sol = adaptive_birkhoff_extrapolation(h, F, x0; Kinit, Kmax, Kstride, Nfactor=5)
        
        # Check rotation vector
        get_w0!(sol, dim)
        @test norm(sol.w0[2:2+dim-1]-w0s[1:dim]) < 1e-10

        # Check parameterization
        adaptive_get_torus!(sol)
        if dim == 1
            @test norm(get_circle_residual((y) -> h(F(hinv(y))), sol.tor, Ntheta)) < 1e-10
        else
            theta_vec = (0:Ntheta-1) .* (2π/Ntheta)
            theta_vecs = [theta_vec for ii = 1:dim]
            @test kam_residual_rnorm(sol.tor, (y) -> h(F(hinv(y))), theta_vecs) < 1e-10
        end
    end


end
