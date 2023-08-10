include("./PeriodicOrbits.jl")
include("./InvariantCircles.jl")
include("./ConnectingOrbit.jl")
include("./Island.jl")
include("./PP_Surrogates.jl")
using Plots
using LinearAlgebra
using Peaks

function s_to_θ(dzs::AbstractArray, θ::AbstractArray, s::AbstractArray)
    i_s = 1;
    i_θ = 1;
    θ_s = zeros(length(s));
    θ_s[end] = θ[end];

    Ns = length(s);
    Nθ = length(dzs);

    s_max = -1+dzs[1];

    while (i_s <= Ns) && (i_θ < Nθ)
        if s[i_s] <= s_max
            θ_s[i_s] = θ[i_θ+1] + (θ[i_θ+1] - θ[i_θ]) * (s[i_s] - s_max)/dzs[i_θ]
            i_s = i_s + 1;
        else
            i_θ = i_θ + 1
            s_max = s_max + dzs[i_θ];
        end
    end

    return θ_s
end

function initialize_connecting_orbit(island::Island, θ0::Number, θ1::Number,
                                      x0::AbstractArray, x1::AbstractArray,
                                      FJ::Function, Na::Integer, p::Integer,
                                      Ns::Integer, is_positively_oriented::Bool)
    #
    z = get_circle(island, get_N_circle(island));

    # Get data for a (probably bad) linear spline approximation of z(θ)
    θ = LinRange(θ0, θ1 > θ0 ? θ1 : 2π+θ1, Ns);
    dzs = [norm(deval(z, θ[i_s]; i_circle=1)) for i_s=1:Ns];
    s_max = sum(dzs[1:end-1]);
    dzs = dzs .*(2/s_max);

    # Sample the outermost circle at chebyshev points
    c = ConnectingOrbit(Na, p);
    s = quad_nodes(c, Ns);
    θ_s = s_to_θ(dzs, θ, s);
    meanθ = (θ_s[end] + θ_s[1])/2;
    θ_s = meanθ .+ 0.99 .* (θ_s .- meanθ)
    z_s = z(θ_s, 1);

    # Scale z_s so that the endpoints are correct
    if x1 != x0
        x_difference = x1 - x0
        z_difference = z_s[:, end] - z_s[:, 1];

        z_s = z_s .* norm(x_difference) / norm(z_difference);

        x_difference = x_difference ./ norm(x_difference)
        z_difference = z_difference ./ norm(z_difference)


        cx = x_difference[1];
        sx = x_difference[2];
        Rx = [cx -sx; sx cx];

        cz = z_difference[1];
        sz = z_difference[2];
        Rz = [cz -sz; sz cz];

        z_s = (Rx*Rz')*z_s
        z_s[1, :] = z_s[1, :] .+ (x0[1] - z_s[1, 1]);
        z_s[2, :] = z_s[2, :] .+ (x0[2] - z_s[2, 1]);
    end

    if !is_positively_oriented
        z_s[:, :] = z_s[:, end:-1:1]
    end


    # Least squares solve for the Chebyshev coefficients
    Φs = get_Φ(c, s);
    set_a!(c, kron(Φs, Matrix(1.0*I, 2, 2))\z_s[:]; i_p=1);

    return c;
end

# Looking for point ii where numerator * ii = 1 mod denominator
function reshuffle_orbit(orbit, numerator, denominator, p)
    if (length(orbit)÷p == 2)
        return copy(orbit)
    end

    ii = 0;
    for jj = 0:denominator-1
        if mod(jj * numerator - 1, denominator) == 0
            ii = 1 + jj*p
        end
    end
    new_orbit = copy(orbit);
    p = size(orbit, 2);
    new_orbit[:, 1:p-ii+1] = orbit[:, ii:p]
    new_orbit[:, p-ii+2:p] = orbit[:, 1:ii-1]

    return new_orbit;
end


function get_orientation(z, FJ, θ0, x0, denominator)
    p = get_p(z);
    A = area(z);

    FJq = FJ_qfold(FJ, p*denominator);
    J = FJq(x0)[2];

    nearest_pt = z(θ0, 1);
    eigenJ = eigen(J);

    P = eigenJ.vectors;
    eigen_loc = P \ (nearest_pt - x0);

    P = P * Diagonal(vec(sign.(eigen_loc)));
    Λ = eigenJ.values

    if abs(Λ[1]) < abs(Λ[2])
        P = P * [0 1; 1 0];
    end


    is_positively_oriented = A * det(P) > 0 ? true : false




    return is_positively_oriented
end

function find_rational_τ(τ_z::Number, denom_max::Integer)
    iota_z = mod(τ_z/(2π), 1);

    iota = 0;
    err = minimum(abs.([1-iota_z, iota_z]));
    numerator = 0;
    denominator = 1;

    for denom = 1:denom_max
        for numer = 1:denom-1
            if gcd(numer, denom) == 1
                iota_i = numer/denom;
                err_i = abs(iota_z - iota_i);
                if err_i < err
                    iota = iota_i;
                    err = err_i;
                    numerator = numer;
                    denominator = denom;
                end
            end
        end
    end

    τ = 2π*iota;

    n_dif = round(Int64, (τ_z - τ)/(2π));
    τ = τ + n_dif*2π;
    numerator = numerator + n_dif * denominator;

    return τ, numerator, denominator
end

function orbit_is_stable(unstable_orbit, FJ)
    d = size(unstable_orbit, 2);
    J = Matrix{Float64}(I, 2, 2);
    for ii = 1:d
        F, Ji = FJ(unstable_orbit[:, ii]);
        J = Ji * J;
    end
    return 4*det(J) > tr(J)^2
end


function find_island_boundary!(island::Island, FJ::Function, Na::Integer, rtol::Number; verbose=false)
    # These are pretty arbitrary values - engineering at its finest
    Nθ = 200;
    w = 5; # Width of minimization for find peaks

    # Get position of suspected unstable equlibria
    z = get_circle(island, get_N_circle(island));
    θ = (1-w:Nθ+w).*(2π/Nθ);
    dz = [norm(deval(z, θi)) for θi in θ]
    inds = argminima(dz, w);

    # Get nearest rational
    τ, numerator, denominator = find_rational_τ(get_τ(z), length(inds))
    set_τ_edge!(island, τ)

    if verbose
        println("Found $(length(inds)) minima. The rotation number is $(numerator)/$(denominator). τ=$(τ)")
    end

    # Find unstable equilibria
    p = get_p(island);
    accepted_inds = copy(inds);
    for ii = 1:length(inds)
        unstable_orbit = zeros(2, p*denominator)
        o_point = get_stable_orbit(island);
        try
            unstable_orbit = newton_periodic(FJ, z(θ[inds[ii]], 1), p*denominator; verbose)
        catch
            for ii = 1:p*denominator
                unstable_orbit[:, ii] = o_point[:, 1]
            end
        end

        # Check if the orbit found is valid
        is_new_orbit = true;
        prev_orbits = get_unstable_orbits(island)
        unstable_orbit_tol = 1e-5 * norm(z);
        for kk = 1:p*denominator
            # If we found the o-point, we are not happy
            if norm(unstable_orbit[:, kk] - o_point[:, 1]) < unstable_orbit_tol
                is_new_orbit = false
            end

            # If we have already found this orbit, we are also not happy
            for jj = 1:length(prev_orbits)
                prev_orbit = prev_orbits[jj];
                if norm(unstable_orbit[:, 1] - prev_orbit[:, kk]) < unstable_orbit_tol
                    is_new_orbit = false
                    orbit_num = jj;
                    orbit_ind = kk;
                end
            end

            if orbit_is_stable(unstable_orbit, FJ)
                is_new_orbit = false;
            end
        end

        # Push the orbit, assuming we are happy with it
        if is_new_orbit
            push_unstable_orbit!(island, unstable_orbit);
        else
            accepted_inds[ii] = -1;
        end

        if verbose
            println("is_new_orbit = $(is_new_orbit), unstable_orbit[ii] = $(unstable_orbit)");
        end
    end

    # Find orientation of stable and unstable directions
    accepted_inds = accepted_inds[accepted_inds .!= -1]
    if length(accepted_inds) == 0
        return
    end

    orbit_1 = get_unstable_orbit(island, 1);
    is_positively_oriented = get_orientation(z, FJ, θ[accepted_inds[1]], orbit_1[:, 1], denominator)
    orientation_found = true
    if verbose
        println("is_positively_oriented = $(is_positively_oriented)")
    end

    # Find connecting orbits
    N_unstable_orbit = get_N_unstable_orbit(island)
    for ii = 1:N_unstable_orbit
        ind0 = accepted_inds[ii]
        ind1 = 1;
        orbit_1 = get_unstable_orbit(island, ii);
        if ii == N_unstable_orbit
            orbit_2 = reshuffle_orbit(get_unstable_orbit(island, 1), numerator, denominator, p)

            err = Inf;
            # Get the closest minimum on z to the next point
            for jj = 1:length(inds)
                err_jj = norm(z(θ[inds[jj]],1) - orbit_2[:, 1]);
                if err_jj < err
                    ind1 = inds[jj];
                    err = err_jj
                end
                if verbose
                    println("x-point initial guess: jj=$jj, ind1=$ind1, err_jj=$err_jj, err=$err, inds[jj]=$(inds[jj]), θ[inds[jj]]=$(θ[inds[jj]])")
                end
            end
        else
            orbit_2 = get_unstable_orbit(island, ii+1);
            ind1 = accepted_inds[ii+1]
        end

        Ns = 2*(Na+1);
        c = initialize_connecting_orbit(island, θ[ind0], θ[ind1], orbit_1[:, 1],
                                        orbit_2[:, 1], FJ, Na, denominator*p, Ns,
                                        is_positively_oriented)
        ends = is_positively_oriented ? hcat(orbit_1[:, 1], orbit_2[:, 1]) : hcat(orbit_2[:, 1], orbit_1[:, 1]);
        gn_connecting!(c, FJ, ends; Ns, verbose, rtol)
        push_connecting_orbit!(island, c)
    end
end


function island_continuation_predict(island::Island, FJ::Function, h_z::Number,
                                     h_τ::Number, h_s::Number, Na::Integer,
                                     max_radius::Number, circle_type::String)
    # Island stuff
    p = get_p(island);
    stable_orbit = get_stable_orbit(island);
    N_circle = get_N_circle(island)

    z = FourierCircle(Na; p);
    c = fixed_θ0_constraint(z);

    if N_circle == 0
        circle_linear!(z, FJ, stable_orbit, h_s*h_z);
        set_τ_center!(island, get_τ(z));
    elseif N_circle == 1
        z_prev = get_circle(island, 1);
        if norm(z_prev) > max_radius
            return 0, 0
        end

        # Δznorm = norm(z_prev);
        Δroot_area = sqrt(abs(area(z_prev))); # Continue in the sqrt of the area
        Δτ = get_τ(z_prev) - get_τ_center(island);
        Δs = sqrt(Δroot_area^2 / h_z^2 + Δτ^2 / h_τ^2);

        set_τ!(z, get_τ(z_prev) + h_s * (Δτ / Δs) )
        if circle_type == "Fourier"
            for i_circle = 1:p
                Δa0 = z_prev[0, i_circle] - stable_orbit[:, i_circle]
                z[0, i_circle] = z_prev[0, i_circle] + h_s .* (Δa0./Δs)
                Na = get_Na(z);
                for i_A = 1:Na
                    ΔAi = z_prev[i_A, i_circle];
                    z[i_A, i_circle] = z_prev[i_A, i_circle] + h_s * (ΔAi./Δs);
                end
            end
        end
    else
        z_prev = get_circle(island, N_circle);
        z_prev2 = get_circle(island, N_circle-1);
        if norm(z_prev) > max_radius
            return 0, 0
        end

        # Δznorm = norm(z_prev) - norm(z_prev2);
        root_area_prev = sqrt(abs(area(z_prev)))
        Δroot_area = root_area_prev - sqrt(abs(area(z_prev2)));
        Δτ = get_τ(z_prev) - get_τ(z_prev2);
        Δa = get_a(z_prev) - get_a(z_prev2)

        Δs = sqrt(Δroot_area^2 / h_z^2 + Δτ^2 / h_τ^2);

        set_τ!(z, get_τ(z_prev) + h_s * (Δτ / Δs) )
        set_a!(z, get_a(z_prev) + h_s .* (Δa ./ Δs))
        root_area_next = sqrt(abs(area(z)));

        # c_A = area_constraint(z)
        # c_τ = fixed_τ_constraint(z)
        # c[2, :] = ((h_z^(-2) * (1 - root_area_prev/root_area_next)) .* c_A +
        #            2*(h_τ^(-2) * (get_τ(z) - get_τ(z_prev))) .* c_τ)
        # c[2, :] = c_A
        # display(c[2,:])
    end


    # constraint = fixed_θ0_constraint; # will be changed
    constraint = (z) -> c;
    return z, constraint
end

function island_continuation_from_center(stable_orbit::AbstractArray,
                                         FJ::Function, h_z::Number, Nθ::Integer,
                                         Na::Integer, max_radius::Number;
                                         maxiter::Integer=10, rtol::Number=1e-8,
                                         constraint::Function=fixed_θ0_constraint,
                                         verbose=false, circle_type="Fourier",
                                         h_τ_min::Number = 0.01, N_refinement::Integer=4,
                                         get_boundary::Bool=true)
    # Get linear island
    island = Island(stable_orbit);

    # Get larger islands until Newton's method doesn't converge or we reach our maximum size


    if verbose
        println("Starting island_continuation_from_center
    h_z = $(h_z), Nθ = $(Nθ), Na = $(Na), max_radius=$(max_radius)
    p = $(get_p(island)), stable_orbit = $(stable_orbit)")
    end

    i_h = 0;
    h_s = 1;
    i_refinement = 0;
    refinement_factor = 2.0;
    h_τ = h_τ_min * refinement_factor^N_refinement;
    # predict
    z, constraint = island_continuation_predict(island, FJ, h_z, h_τ, h_s, Na,
                                                max_radius, circle_type)
    while z != 0
        # correct
        niter = gn_circle(FJ, z; Nθ, maxiter, rtol, verbose, λ=0.0, constraint)

        if niter > maxiter
            if i_refinement < N_refinement
                # refine
                i_refinement = i_refinement + 1;
                h_s = h_s / refinement_factor;
            else
                # give up
                if get_boundary
                    find_island_boundary!(island, FJ, Na, rtol; verbose)
                end

                if verbose
                    println("Exiting continuation")
                end
                return island
            end

            if verbose
                i_h = i_h+1;
                println("i_h = $(i_h), i_refinement = $(i_refinement)")
            end
        else
            if (niter > (2*maxiter)÷3) && (i_refinement < N_refinement)
                # refine if convergence took too long
                i_refinement = i_refinement + 1;
                h_s = h_s / refinement_factor;
            end
            push_circle!(island, z);

            if verbose
                i_h = i_h+1;
                println("i_h = $(i_h), i_refinement = $(i_refinement), radius = $(norm(z)), τ = $(get_τ(z))")
            end
        end

        # predict
        z, constraint = island_continuation_predict(island, FJ, h_z, h_τ, h_s,
                                                    Na, max_radius, circle_type)
    end

    return island;
end
