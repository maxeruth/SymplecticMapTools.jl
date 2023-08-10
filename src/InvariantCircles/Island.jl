# Collection of invariant circles comprising an island
# Stored in an doubling array
mutable struct Island
    stable_orbit::AbstractArray;      # Stable point at the center of the island
    τ_center::Number;                 # Rotation number at magnetic axis
    τ_edge::Number;                   # Rotation number at the island edge

    circles::AbstractArray;           # Invariant circles
    N_circle::Integer;                # Number of circles currently stored

    unstable_orbits::AbstractArray;   # unstable points at boundary of circle
    N_unstable_orbit::Integer;        # Number of unstable orbits

    connecting_orbits::AbstractArray; # segmants connecting the unstable points
    N_connecting_orbit::Integer;      # Number of connecting orbits

    p::Integer;                       # Period

    function Island(stable_orbit::AbstractArray)
        p = length(stable_orbit)÷2;
        new(reshape(deepcopy(stable_orbit), 2, p), 0.0, 0.0, [], 0, [], 0, [], 0, p);
    end
end

function get_p(island::Island)
    return island.p;
end

function get_stable_orbit(island::Island)
    return island.stable_orbit[:,:];
end

function get_τ_center(island::Island)
    return island.τ_center;
end

function set_τ_center!(island::Island, τ_center::Number)
    island.τ_center = τ_center;
end

function get_τ_edge(island::Island)

    return island.τ_edge;
end

function set_τ_edge!(island::Island, τ_edge::Number)
    island.τ_edge = τ_edge;
end


function get_N_circle(island::Island)
    return island.N_circle;
end

function set_N_circle!(island::Island, N_circle::Integer)
    island.N_circle = N_circle;
end

function get_circles(island::Island)
    return island.circles[1:get_N_circle(island)];
end

function get_circle(island::Island, ii::Integer)
    return island.circles[ii];
end

function set_circles!(island::Island, circles::AbstractArray)
    island.circles = circles;
end

function push_circle!(island::Island, circle::InvariantCircle)
    N_circle = get_N_circle(island)

    # If necessary, double the size of the array
    if N_circle == length(island.circles)
        doubled_array = Vector{InvariantCircle}(undef, N_circle == 0 ? 1 : 2*N_circle);

        doubled_array[1:N_circle] = get_circles(island);
        set_circles!(island, doubled_array);
    end

    # put in the new circle, update the number of circles
    island.circles[N_circle+1] = circle;
    set_N_circle!(island, N_circle+1);
end

function get_N_unstable_orbit(island::Island)
    return island.N_unstable_orbit;
end

function set_N_unstable_orbit!(island::Island, N_unstable_orbit)
    island.N_unstable_orbit = N_unstable_orbit;
end

function get_unstable_orbits(island::Island)
    return island.unstable_orbits[1:get_N_unstable_orbit(island)];
end

function get_unstable_orbit(island::Island, i_orbit::Integer)
    return island.unstable_orbits[i_orbit];
end

function set_unstable_orbits!(island::Island, unstable_orbits::AbstractArray)
    island.unstable_orbits = unstable_orbits;
end

function push_unstable_orbit!(island::Island, unstable_orbit::AbstractArray)
    N = get_N_unstable_orbit(island)

    # If necessary, double the size of the array
    if N == length(island.unstable_orbits)
        doubled_array = Vector{Array}(undef, N == 0 ? 1 : 2*N);

        doubled_array[1:N] = get_unstable_orbits(island);
        set_unstable_orbits!(island, doubled_array);
    end

    island.unstable_orbits[N+1] = unstable_orbit;
    set_N_unstable_orbit!(island, N+1);
end

function get_N_connecting_orbit(island::Island)
    return island.N_connecting_orbit;
end

function set_N_connecting_orbit!(island::Island, N_connecting_orbit::Integer)
    island.N_connecting_orbit = N_connecting_orbit
end

function get_connecting_orbits(island::Island)
    return island.connecting_orbits[1:get_N_connecting_orbit(island)]
end

function get_connecting_orbit(island::Island, ii::Integer)
    return island.connecting_orbits[ii];
end

function set_connecting_orbits!(island::Island, connecting_orbits::AbstractArray)
    island.connecting_orbits = connecting_orbits;
end

function push_connecting_orbit!(island::Island, connecting_orbit::ConnectingOrbit)
    N = get_N_connecting_orbit(island)

    # If necessary, double the size of the array
    if N == length(island.connecting_orbits)
        doubled_array = Vector{ConnectingOrbit}(undef, N == 0 ? 1 : 2*N);

        doubled_array[1:N] = get_connecting_orbits(island);
        set_connecting_orbits!(island, doubled_array);
    end

    island.connecting_orbits[N+1] = connecting_orbit;
    set_N_connecting_orbit!(island, N+1);
end


function area(island::Island; i_p::Integer=1, rtol = 1e-8)
    p = get_p(island);
    stable_orbit = get_stable_orbit(island);
    origin = stable_orbit[:, 1];
    A = 0.0;
    for c = get_connecting_orbits(island)
        for ii = i_p:p:get_p(c)
            A = A + area(c; origin, i_p=ii, rtol)
        end
    end

    return A
end



function Plots.plot!(p::Plots.Plot, island::Island; d_circle::Integer=1, N=50, circ_linewidth=1, conn_linewidth=2)
    for ii = d_circle:d_circle:get_N_circle(island)
        z = get_circle(island, ii)
        Plots.plot!(p, z; linewidth=circ_linewidth, color=:black, N);
    end

    stable_orbit = get_stable_orbit(island);
    Plots.scatter!(stable_orbit[1, :], stable_orbit[2, :], color=:blue, label=false, markersize=7)


    for c in get_connecting_orbits(island)
        Plots.plot!(c, color=:purple; linewidth=conn_linewidth, N)
    end

    for orbit in get_unstable_orbits(island)
        Plots.scatter!(orbit[1, :], orbit[2, :], marker=:x, color=:blue, label=false, markersize=7)
    end
end

function Plots.plot(island::Island; d_circle::Integer=1, N=50, circ_linewidth=1, conn_linewidth=2)
    p = Plots.plot();
    Plots.plot!(p, island; d_circle, N, circ_linewidth, conn_linewidth);

    return p;
end

function plot_area_vs_iota(island::Island; xlabel="Area", ylabel="Rotational Transform (ι)", is_boundary=true, root_area=false)
    N_circle = get_N_circle(island);
    N_connecting_orbit = get_N_connecting_orbit(island);
    is_boundary = is_boundary && (N_connecting_orbit > 0);
    areas = zeros(is_boundary ? 2+N_circle : 1+N_circle);
    iotas = zeros(is_boundary ? 2+N_circle : 1+N_circle);

    iotas[1] = get_τ_center(island)/(2π);

    for i_circle = 1:N_circle
        z = get_circle(island, i_circle);
        areas[i_circle+1] = abs(area(z));
        iotas[i_circle+1] = get_τ(z)/(2π);
    end

    if is_boundary
        areas[end] = abs(area(island));
        iotas[end] = get_τ_edge(island) / (2π)
    end

    if root_area
        areas[:] = sqrt.(areas)
        if xlabel == "Area"
            xlabel = "Root Area"
        end
    end

    p = Plots.scatter(areas, iotas, label=false);
    Plots.xlabel!(xlabel);
    Plots.ylabel!(ylabel);

    return p;
end
