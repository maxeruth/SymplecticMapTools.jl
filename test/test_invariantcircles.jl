@testset "InvariantCircles.jl and FourierCircle.jl" begin
    k_sm = 0.2
    FJ_sm = standard_map_FJ(k_sm)
    h, HJ, hinv, HJinv = polar_map(;z0 = -0.2)
    FJ = (y) -> begin
        x, Jx = HJinv(y)
        F, J_sm = FJ_sm(x)
        J = J_sm*Jx
        y_next, Jy_next = HJ(F)
        (y_next, Jy_next*J)
    end

    y0 = zeros(2,2);
    y0[:, 1] = h([0.5, 0.5]);
    y0[:, 2] = h([1.0, 0.5])
    # Test that the orbit is periodic
    @test norm(y0[:,1] - FJ(y0[:,2])[1]) < 1e-12
    @test norm(y0[:,2] - FJ(y0[:,1])[1]) < 1e-12

    # Finite difference check the derivative
    dy = [0.1, 0.2]*1e-4;
    F, J = FJ(y0[:,1]+dy);
    @test norm((y0[:, 2] + J*dy) - F) < 1e-7


    # Test getter/setters for FourierCircle
    Na = 10
    p = 2
    z =  FourierCircle(Na; p)
    @test get_Na(z) == Na
    @test get_p(z) == p
    z[0, 1] = [1., 1.]
    @test z[0,1] == [1., 1.]
    z[1, 2] = [1., 0., 0., 1.]
    @test z[1,2] == [1., 0., 0., 1.]
    set_τ!(z, π)
    @test get_τ(z)-π == 0.

    # Test evaluation stuff
    @test norm(z(0., 2) - [1., 0.]) < 1e-14
    @test norm(deval(z, 0.; i_circle=2)-[0., 1.]) < 1e-14
    zp = deriv(z)
    @test norm(zp(0., 2)-[0., 1.]) < 1e-14
    @test abs(area(z; i_circle=2) - π) < 1e-14
    @test norm(shifted_eval(z, [π]; i_circle=2) - [1., 0.]) < 1e-14

    # Get island around the orbit
    h_circ = 0.01
    z =  FourierCircle(Na; p)
    circle_linear!(z, FJ, y0, h_circ)
    @test abs(average_radius(z, i_circle=1)-h_circ) < 1e-12

    # Newton iterate the island
    Nθ = 100;
    maxiter = 10;
    n = gn_circle(FJ, z, Nθ; maxiter)
    @test n < maxiter

    # Check residual
    F_sm = standard_map_F(k_sm)
    F = (y) -> h(F_sm(hinv(y)))
    res = get_circle_residual(F, z, Nθ)
    @test norm(res) < 1e-10
end
