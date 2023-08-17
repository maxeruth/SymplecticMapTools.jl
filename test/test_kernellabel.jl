
@testset "Kernel Label" begin
    # Create Kernel Label
    kernels = [:SquaredExponential, :InverseMultiquadric, :FourierSE, :Fourier]

    x = [0.0; 0.0;; 0.0; 0.5;; 0.5; 0.5;; 0.5; 0.0]
    c = [1.0, 0.0];
    σ = 0.5;

    # Test positive definiteness on 2x2 matrix for every kernel
    for (ii, kernel) in enumerate(kernels)
        k = KernelLabel(x, c, σ; kernel)
        A = get_matrix(k, x)[1:2, 1:2]
        @test norm(A'-A) ≈ 0.
        @test tr(A) > 0
        @test det(A) > 0
    end

    # Test windowing functions
    lims = [-0.25, 0.25]
    α = 1e-10
    w = window_weight(x, lims, α)
    @test norm(w - [0.0, 1.0, 1.0, 0.0]) ≈ 0.
    w_rec = rectangular_window_weight(x, lims, lims, α)
    @test norm(w_rec - [0.0, 1.0, 1.0, 1.0]) ≈ 0.

    # Sample the standard map
    k_sm = 0.1;
    N = 100;
    xb = [0., 1.]
    yb = [0., 1.]
    F = standard_map_F(k_sm);

    xs = kernel_sample_F(F, N, xb, yb)
    @test norm(F(xs[:, 1])-xs[:,2]) ≈ 0.0

    # Test kernel_eigs
    ϵ = 1e-8;
    nev = 5;
    σ = 0.5;
    kernel = :FourierSE;
    lims = [0.05, 0.95];
    α = 1e-4
    w = window_weight(xs, lims, α)

    λs, vs, k_eigs = kernel_eigs(xs, ϵ, nev, σ, w; kernel)
    @test abs(λs[1]/4.85624090981917e-7 - 1) < 1e-4 # This is arbitrary, but it makes sure the behavior doesn't change at least...

    # And kernel_bvp
    boundary_values = [tanh((2xi[2]-1)/α) for xi = eachcol(xs)]
    k, R, R_bc, R_inv, R_eps = kernel_bvp(xs, ϵ, σ, w, boundary_values; kernel)
    @test (R/3.430457448227529e-6 - 1) < 1e-4 # Again, arbitrary

    # Test that get_energies is returning the same stuff as kernel_bvp
    EK, EInv, Ebd, EL2 = get_energies(k; W = Diagonal(w))
    @test abs(R_inv - EInv) < 1e-4
    @test (R_eps/ϵ - EK) < 1e-4
end
