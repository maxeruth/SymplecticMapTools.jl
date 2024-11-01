@testset "Lyapunov Exponents" begin
    x0 = [0.0,0.1]
    k_sm = 0.2
    FJ_sm = standard_map_FJ(k_sm)
    N = 1000
    L_vec = lyapunov_exponents(x0,FJ_sm,N)
    L = lyapunov_exponent(x0,FJ_sm,N)

    @test L == L_vec[end]
    @test abs(L) < 0.01
end