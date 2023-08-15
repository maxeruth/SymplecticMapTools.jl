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
    c, sums, resid, xs, hs, rnorm, K, history = adaptive_birkhoff_extrapolation(
       h, F, x0; Kinit, Kmax, Kstride, iterative=true, Nfactor=1)
    cls, sumsls, residls, xsls, hsls, rnormls, Kls, historyls = adaptive_birkhoff_extrapolation(
       h, F, x0; Kinit, Kmax, Kstride, iterative=false, Nfactor=1.1)

    # Test the iterative and direct methods give similar predictions
    @test norm(sums[:,1] - sumsls[:,1]) < 1e-8
    @test norm(xs[:, 1:Kinit] - xsls[:, 1:Kinit]) < 1e-14
    @test norm(hs[:, 1:Kinit] - hsls[:, 1:Kinit]) < 1e-14

    # Birkhoff average should match the rotation number
    @test abs(dot(xs[2,1:2K-1],c) - 0.5) < 1e-8

    # Test the invariant circle is correct
    Nθ = 100

    z = get_circle_info(hs, c)
    @test get_p(z) == 2
    @test norm(get_circle_residual((y) -> h(F(hinv(y))), z, Nθ)) < 1e-10
end
