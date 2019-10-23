using Oceananigans, PyPlot, Random, Printf

using Oceananigans: Face, Cell

# Resolution
Nx = 64
Ny = 64
Nz = 64 

# Domain size
Lx = 1000e3
Ly = 1000e3
Lz = 800

# Physical parameters
 f = 1e-4     # Coriolis parameter
N² = 1e-8     # Stratification in the "halocline"
 h = 200      # Halocline depth / extent
Δb = N² * h

const τ₀ = -1e-5    # Wind stress magnitude
const Δτ = 80e3     # Width of wind stress

# Simulation end time
end_time = 30day

# Wind stress / velocity flux forcing functions
@inline τˣ(x, y, t) = τ₀ * y/Δτ * exp(-(x^2 + y^2) / 2Δτ^2)
@inline τʸ(x, y, t) = τ₀ * x/Δτ * exp(-(x^2 + y^2) / 2Δτ^2)

u_bcs = HorizontallyPeriodicBCs(   top = FunctionBoundaryCondition(Flux, :z, Face, Cell, τˣ),
                                bottom = BoundaryCondition(Value, 0))

v_bcs = HorizontallyPeriodicBCs(   top = FunctionBoundaryCondition(Flux, :z, Cell, Face, τʸ),
                                bottom = BoundaryCondition(Value, 0))

b_bcs = HorizontallyPeriodicBCs(   top = BoundaryCondition(Value, Δb),
                                bottom = BoundaryCondition(Value, 0))

model = Model(
         architecture = CPU(),
                 grid = RegularCartesianGrid(size=(Nx, Ny, Nz), x=(-Lx/2, Lx/2), y=(-Ly/2, Ly/2), z=(-Lz, 0)),
             coriolis = FPlane(f=f),
             buoyancy = BuoyancyTracer(),
              tracers = :b,
              closure = ConstantAnisotropicDiffusivity(νv=1e-1, νh=1, κv=1e-1, κh=1),
  boundary_conditions = BoundaryConditions(u=u_bcs, v=v_bcs, b=b_bcs)
)

# Initial condition
b₀(x, y, z) = Δb * exp(z^2 / 2h^2)
set!(model, b=b₀)

wizard = TimeStepWizard(cfl=0.005, Δt=minute, max_change=1.1, max_Δt=5minute)

# A diagnostic that returns the maximum absolute value of `w` by calling
# `wmax(model)`:

wmax = FieldMaximum(abs, model.velocities.w)

# Set up output
fields_to_output = merge(model.velocities, (b=model.tracers.b,))
output_writer = JLD2OutputWriter(model, FieldOutputs(fields_to_output); interval=day, prefix="beaufort_gyre",
                                 force=true)

# Create a figure
fig, axs = subplots(ncols=2, figsize=(12, 6))

function makeplot!(fig, axs, model)

    fig.suptitle("\$ t = \$ $(prettytime(model.clock.time))")

    sca(axs[1]); cla()
    title("\$ u(x, y, z=0) \$")
    imshow(interior(model.velocities.u)[:, :, Nz])

    kplot = Nz - 2
    sca(axs[2]); cla()
    title("\$ w(x, y, z=$(model.grid.zF[kplot]))")
    imshow(interior(model.velocities.w)[:, :, kplot])
    
    [ax.tick_params(left=false, labelleft=false, bottom=false, labelbottom=false) for ax in axs]

    return nothing
end

# Finally, we run the the model in a `while` loop.

while model.clock.time < end_time

    ## Update the time step associated with `wizard`.
    update_Δt!(wizard, model)

    ## Time step the model forward
    walltime = @elapsed time_step!(model, 10, wizard.Δt)

    ## Print a progress message
    @printf("i: %04d, t: %s, Δt: %s, wmax = %.1e ms⁻¹, wall time: %s\n",
            model.clock.iteration, prettytime(model.clock.time), prettytime(wizard.Δt),
            wmax(model), prettytime(walltime))

    model.architecture == CPU() && makeplot!(fig, axs, model)
end
