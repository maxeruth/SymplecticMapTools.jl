using SymplecticMapTools   # Used for finding invariant tori
using GeometricIntegrators # Used for symplectic integration
using CairoMakie           # Used for plotting

# # An example of 3D invariant tori
# Here, we are showing an example of three three-dimensional invariant tori computed via the Birkhoff RRE process.
# Through this, we reproduce a figure from the [Ruth, Kulik, and Burby](https://doi.org/10.48550/arXiv.2505.08715).
# The system we consider is the elliptic three-body problem, described by the Hamiltonian
# ```math
# H(\mathbf x, \mathbf p, f) = \frac{1}{2}\left(\Vert\mathbf x\Vert^2+\Vert\mathbf p\Vert^2\right) + y p_x - x p_y - \phi(\mathbf x,f), 
# ```
# where $\mathbf x = (x,y,z)$ is the position of a light object (such as a satellite), $\mathbf p = (p_x, p_y, p_z)$ is the momentum in a co-rotating elliptic frame, $f$ is the true anomaly, and the potential $\phi$ is
# ```math
#     \phi = \frac{1}{1+ \epsilon \cos f}\left[\frac{1}{2}\left((1-\mu)\rho_1^2 + \mu \rho_2^2\right) + \frac{1-\mu}{\rho_1} + \frac{\mu}{\rho_2} \right].
# ```
# The potential $\phi$ is the gravitational pull of the two primary bodies with mass ratio $\mu$ and eccentricity $\epsilon$ chosen to match the earth-moon system, with the $\rho_j$ factors defined by
# ```math
# \rho_1^2 = (x-\mu)^2 + y^2+z^2, \quad \rho_2^2 = (x - (\mu-1))^2 + y^2+z^2.
# ```
# We obtain a symplectic map by symplectically integrating the Hamiltonian dynamical system
# ```math
# \frac{d\mathbf x}{df} = -\frac{\partial H}{\partial \mathbf p}, \qquad \frac{d\mathbf p}{df} = \frac{\partial H}{\partial \mathbf x},
# ```
# over an orbit of the moon about the earth, i.e. as $f$ increases from $0$ to $2\pi$.
# 
# To begin, we define the path where we will save our results and the dynamical system.

## Location to save output to (saves to path/of/SymplecticMapTools.jl/docs/3D_ER3BP/ by default)
path = dirname(pathof(SymplecticMapTools))*"/../docs/3D_ER3BP/"
if !isdir(path)
    mkdir(path)
end

#-

## Potential energy due to a single primary
function get_ρ(x,y,z,μ)
    sqrt((x-μ)^2 + y^2 + z^2)
end

## Value of potential energy w.r.t. coordinates
function gradρ(x,y,z,μ)
    dx = x-μ
    ρ = sqrt(dx^2 + y^2 + z^2)
    dρdx = dx/ρ
    dρdy = y/ρ
    dρdz = z/ρ
    return ρ,dρdx,dρdy,dρdz
end

## ER3BP potential gradient w.r.t. ρi
function gradω(ρ1,ρ2,f,ecc,μ)
    dΩdρ1 = (1. - μ)*(ρ1 - 1. / ρ1^2)
    dΩdρ2 = μ*(ρ2 - 1. / ρ2^2)
    prefactor = 1. / (1. + ecc*cos(f))

    return prefactor*dΩdρ1, prefactor*dΩdρ2
end

## ER3BP potential energy
function get_ω(ρ1,ρ2,f,ecc,μ)
    (1. / (1. + ecc*cos(f)) ) * ( 0.5 * ((1-μ)*ρ1^2 + μ*ρ2^2) + (1-μ)/ρ1 + μ/ρ2) 
end

## ER3BP velocity
function v_ER3BP(v, f, q, p, params)
    x,  y, z = q
   px, py, pz = p
   
   ecc = params[:ecc]
   μ   = params[:μ]
   
   v[1] = px+y
   v[2] = py-x
   v[3] = pz
end

## ER3BP force
function f_ER3BP(v, f, q, p, params)
    x,  y, z = q
   px, py, pz = p
   
   ecc = params[:ecc]
   μ   = params[:μ]

   
   ρ1,dρ1dx,dρ1dy,dρ1dz=gradρ(x,y,z,-μ)
   ρ2,dρ2dx,dρ2dy,dρ2dz=gradρ(x,y,z,1. - μ)
   dωdρ1,dωdρ2=gradω(ρ1,ρ2,f,ecc,μ)

   v[1] = -x + py + dωdρ1*dρ1dx + dωdρ2*dρ2dx
   v[2] = -y - px + dωdρ1*dρ1dy + dωdρ2*dρ2dy
   v[3] = -z + dωdρ1*dρ1dz + dωdρ2*dρ2dz
end

## ER3BP Hamiltonian
function H_ER3BP(f, q, p, params)
    x,  y, z = q
   px, py, pz = p
   ecc = params[:ecc]
   μ   = params[:μ]

   ρ1 = get_ρ(x, y, z,  -μ     )
   ρ2 = get_ρ(x, y, z, 1. - μ)
   ω  = get_ω(ρ1,ρ2,f,ecc,μ)

   0.5*(px^2 + py^2 + pz^2 + x^2 + y^2 + z^2) + y*px - x*py - ω
end


## ER3BP map (parameters default to the earth)
function F_ER3BP(x0::AbstractVector; ecc=.0549, μ=1.215058560962404e-2, Nstep::Integer=1000, integrator=GeometricIntegrators.SRK3())
    @assert length(x0) == 6
    
    ## Convert state vector to (q,p)
    q0 = [x0[1],x0[3],x0[5]]
    p0 = [x0[2]-x0[3],x0[4]+x0[1],x0[6]]

    ## Flow the state
    tspan = (0,2π)
    params = (ecc=ecc, μ=μ)
    prob = HODEProblem(v_ER3BP, f_ER3BP, H_ER3BP, tspan, (1-eps())*(2π)/Nstep, q0, p0; parameters=params)
    int = GeometricIntegrator(prob, integrator)
    sol = integrate(int)
    
    ## Convert back to state vector
    q = sol.q[end]
    p = sol.p[end]
    [q[1], p[1] + q[2], q[2], p[2] - q[1], q[3], p[3]]
end


# Below, we give three different initial conditions which can be used to find invariant tori. 
# The first is a librational orbit near the $L_4$ Lagrange point, the second is a western low prograde orbit, and the third is a distant retrograde orbit.
# These are the orbits displayed in [Ruth, Kulik, and Burby](https://doi.org/10.48550/arXiv.2505.08715).
# We note that each invariant torus will take a while to compute (about 5 minutes on an Apple M2 chip), due primarily to the cost of using a symplectic integrator to integrate the individual trajectories.

## Map and observation function
F = (x) -> F_ER3BP(x)
h = (x) -> x;

## Initial point
x0, file = ([0.5-1.215058560962404e-2+.001, -.002, sqrt(3/4)+.002, .002, .05, .05], "Trojan.jld2") # L4 Librational
## x0, file = ([1.0425396498587201,0,0,.43016606592111284+.001,.001,.001], "Western_Low_Prograde.jld2") # Western Low Prograde
## x0, file = ([9.5460504923700373e-1+.001, +.001, +.001, 6.3987767952607477e-1+.001,+.001,+.001], "Distant_Retrograde.jld2") # Distant retrograde


### Birkhoff RRE parameters
rtol = 1e-13;  # RRE tolerance
Kinit = 500;   # Filter size (reducing Kmax speeds up the algorithm)
Kmax = 2500;
Kstride = 500;
Nfactor = 4;   # rectangularity of RRE least-squares
d = 3          # dimensionality of torus

## Load or compute the invariant torus
sol = nothing
if isfile(path*file)
    println("Loading Solution")
    sol = load_rre(path*file);
else
    ## Compute the RRE Filter
    println("Computing Solution (this will take a few minutes)")
    sol = adaptive_birkhoff_extrapolation(h, F, x0; rtol, Kinit, Kmax, Kstride, Nfactor);
    println("The extrapolation returned a filter with size $(sol.K) and RRE error $(sol.resid_rre)"); flush(stdout)
    
    ## Find the frequency
    get_w0!(sol, d; Nsearch=20)
    println("The rotation vector was found to be $(sol.w0[2:end]) with log posterior $(sol.resid_w0)"); flush(stdout)
    
    ## Find the torus
    adaptive_get_torus!(sol) 
    println("The torus was found with resolution $(size(sol.tor.a)[3:end]) with validation residual $(sol.resid_tor)"); flush(stdout)
    
    ## Find the KAM torus residual (this is a good a posteriori check, but takes some time)
    Ntheta = 11 # Number of points to check in each torus direction
    thetavecs = [(0:Ntheta-1) .* (2π/Ntheta) for ii = 1:d]
    resid_kam = kam_residual_rnorm(sol.tor, F, thetavecs) # Get the residual (this can be slow as well)
    println("The torus KAM residual was $(resid_kam)"); flush(stdout)

    ## Save the result
    save_rre(path*file, sol)
end
;

#-

# Now, we plot the torus we computed.
# Because it is a 3D torus, we plot 2D torus slices.
# If one or both of the other two invariant torus cases are uncommented and run, then we plot those tori as well.

#-

## Get 3D grid of points to plot on the torus
function plot_eval_on_grid(tor, Nθs, Q)
    θvec = [(0:Nθ-1) .* (2π/Nθ) for Nθ in Nθs]
    x = evaluate_on_grid(tor, θvec);
    Nisland = size(x,2)
    
    x_plot = zeros(size(Q,2),Nisland,Nθs...)
    for ii = 1:Nθs[1], jj = 1:Nθs[2], kk = 1:Nθs[3]
        x_plot[:,:,ii,jj,kk] = Q'x[:,:,ii,jj,kk]
    end
    x_plot
end

## Mesh and plot an invariant torus
function torus_mesh!(ax, xs; linewidth=1, color=(:black,0.05), kwargs...)
    Ns = size(xs)[3:4]
    Nisland = size(xs,2)
    D = size(xs,1)
    
    faces = zeros(Integer, 3, 2, Ns[1], Ns[2], Nisland)
    for ii = 1:Nisland
        for jj = 1:Ns[1], kk = 1:Ns[2]
            j1 = jj
            j2 = mod1(jj+1,Ns[1])
            k1 = kk
            k2 = mod1(kk+1,Ns[2])

            l1 = (k1-1) * Ns[1] + j1
            l2 = (k1-1) * Ns[1] + j2
            l3 = (k2-1) * Ns[1] + j1
            l4 = (k2-1) * Ns[1] + j2
            
            faces[:,1,jj,kk,ii] = [l1, l2, l3]
            faces[:,2,jj,kk,ii] = [l4, l3, l2]
        end
        mesh!(reshape(xs[:,ii,:,:], D, Ns[1]*Ns[2]), reshape(faces[:,:,:,:,ii], 3, 2*Ns[1]*Ns[2])', color = color)
    end
end;

#- 

## Check whether we have computed the given torus
files = ["Trojan.jld2", "Western_Low_Prograde.jld2","Distant_Retrograde.jld2"]
orbits = [L"$L_4$ librational", L"$$Western low prograde", L"$$Distand retrograde"]
anglelabels = [L"$\theta_1$ slices",L"$\theta_2$ slices",L"$\theta_3$ slices"]
isfiles = [isfile(path*file) for file in files]
Nfiles = sum(isfiles)


## Projection to plot the figure onto
Q = zeros(6,3) # 3D Projection matrix
dims = [1,3,5] # Dimensions to project onto. (1,3,5) are position (x,y,z) coordinates, (2,4,6) are velocity (v_x,v_y,v_z)
labels = [L"x",L"v_x",L"y",L"v_y",L"z",L"v_z"]
Q[dims[1],1]=1.
Q[dims[2],2]=1.
Q[dims[3],3]=1.


## Plot the moon as a point?
plot_moon = [false,true,true]
mu = 1.215058560962404e-2

## Color gradient and resolution to plot the tori in
Ncolor = 8
cticks = (0:7).*(2π/8)
cticklabels = [L"0",L"\pi/4",L"\pi/2",L"3\pi/4",L"\pi",L"5\pi/4", L"3\pi/2", L"7\pi/4"]
cg = cgrad(:Paired_8, Ncolor, categorical=true)
Nθss = [[Ncolor,50,50],[50,Ncolor,50],[50,50,Ncolor]]

## Make the figure
f = Figure(size=(1200,100 + (950-100)*(Nfiles/3)))
colgap!(f.layout, 0) 
rowgap!(f.layout, 0)
elevations = [[11π/32,π/8,3π/8],[5π/16,31π/64,3π/8],[5π/16,31π/64,3π/8]] # Elevation of the plots
prot = 12
zticks = [([-5e-2,0,5e-2],["-5e-2","0","5e-2"]),([-1e-3,0,1e-3],["-1e-3","0","1e-3"]),([-1e-3,0,1e-3],["-1e-3","0","1e-3"])]

for ifile in (1:3)[isfiles] # Loop over orbit type
    soli = load_rre(path*files[ifile])
    tor = soli.tor
    Label(f[ifile, 0], orbits[ifile], rotation = pi/2,tellheight=false,tellwidth=true)
    for lp = 1:3 # Loop over slicing direction
        if ifile == 1
            Label(f[0, lp], anglelabels[lp],tellheight=true,tellwidth=false) # Label of the orbit type
        end

        ## Create 3D Axis
        ax = Axis3(f[ifile,lp],xlabel=labels[dims[1]],
                   ylabel=labels[dims[2]],zlabel=labels[dims[3]],
                   elevation=elevations[ifile][lp], azimuth=2π/16-π/2, 
                   zticks=zticks[ifile], protrusions = prot)
        

        ## Evaluate the torus
        Nθs = Nθss[lp]
        xp = plot_eval_on_grid(tor, Nθs, Q);

        ## Plot the moon
        if plot_moon[ifile]
            scatter!([1-mu],[0.],[0.];color=:black,markersize=20)
        end

        
        for ii = 1:Nθs[lp] # Loop over slices
            color = (cg[ii],0.2) # add some transparency to the plot
            ## plot the torus slices
            if lp == 1
                torus_mesh!(ax, xp[:,:,ii,:,:]; color=color)
            elseif lp == 2
                torus_mesh!(ax, xp[:,:,:,ii,:]; color=color)
            else
                torus_mesh!(ax, xp[:,:,:,:,ii]; color=color)
            end
        end
    end
end

Colorbar(f[:,4],colormap=cg, limits = (-π/8,2π - (π/8)), ticks=(cticks,cticklabels),
         label=L"θ_j")

save(path*"Tori.png", f)

f