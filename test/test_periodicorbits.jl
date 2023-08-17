@testset "PeriodicOrbits.jl" begin
    k_sm = 0.2
    FJ = standard_map_FJ(k_sm)
    x0 = [0.49, 0.49];
    q = 2;

    xs_BFGS, res = BFGS_periodic(FJ, x0, q)
    xs_newton, res_newton = newton_periodic(FJ, x0, q)
    xs = [0.5; 0.5;; 1.0; 0.5]

    @test norm(xs-xs_BFGS) < 1e-7
    @test norm(xs-xs_newton) < 1e-7
end
