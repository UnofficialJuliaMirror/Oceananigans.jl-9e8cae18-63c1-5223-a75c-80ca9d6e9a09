using Oceananigans, PyPlot, Random, Printf

Nx = 64 
Ny = 64 
Nz = 32 

Lx = 1000e3
Ly = 1000e3
Lz = 800

 f = 1e-4
 α = 0.0      # Thermal expansion coefficient
 β = 8e-4     # Haline contraction coefficient

surface_salinity = 10
 bottom_salinity = 20

const τ₀ = 1e-4     # Kinematic wind stress magnitude
const Δτ = 100e3     # Kinematic wind stress magnitude

τˣ(x, y, t) = τ₀ * y/Δτ * exp(-y^2 / 2Δτ^2)
τʸ(x, y, t) = τ₀ * x/Δτ * exp(-x^2 / 2Δτ^2)

@inline τˣ_kernel(i, j, grid, time, args...) = @inbounds τˣ(grid.xF[i], grid.yC[j], time)
@inline τʸ_kernel(i, j, grid, time, args...) = @inbounds τʸ(grid.xC[i], grid.yF[j], time)

u_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, τˣ_kernel),
                                bottom = BoundaryCondition(Value, 0))

v_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, τʸ_kernel),
                                bottom = BoundaryCondition(Value, 0))

S_bcs = HorizontallyPeriodicBCs(   top = BoundaryCondition(Value, surface_salinity),
                                bottom = BoundaryCondition(Value, bottom_salinity))

model = Model(
         architecture = CPU(),
                 grid = RegularCartesianGrid(size=(Nz, Nz, Nz), length=(Lx, Ly, Lz)),
             coriolis = FPlane(f=f),
             buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=α, β=β)),
              closure = ConstantAnisotropicDiffusivity(νv=1e-4, νh=1e-2, κv=1e-4, κh=1e-2),
  boundary_conditions = BoundaryConditions(u=u_bcs, v=v_bcs, S=S_bcs)
)

## Temperature initial condition: a stable density tradient with random noise superposed.
smooth_step(x, Δ) = (tanh(x/Δ) + 1) / 2
initial_salinity(x, y, z) = bottom_salinity + (surface_salinity - bottom_salinity) * smooth_step(z+100, 20)

set!(model, S=initial_salinity)

wizard = TimeStepWizard(cfl=0.2, Δt=minute, max_change=1.1, max_Δt=5.0)

# A diagnostic that returns the maximum absolute value of `w` by calling
# `wmax(model)`:

wmax = FieldMaximum(abs, model.velocities.w)

# We also create a figure and define a plotting function for live plotting of results.

fig, axs = subplots(ncols=3, figsize=(12, 5))


"""
    makeplot!(axs, model)

Make a triptych of x-z slices of vertical velocity, temperature, and salinity
associated with `model` in `axs`.
"""
function makeplot!(axs, model)
    jhalf = floor(Int, model.grid.Nz/2)

    kinetic_energy_slice(u, v) = @. @views (u[:, jhalf, :]^2 + v[:, jhalf, :]^2) / 2

    ## Coordinate arrays for plotting
    xC = repeat(model.grid.xC, 1, model.grid.Nz)
    zF = repeat(reshape(model.grid.zF[1:end-1], 1, model.grid.Nz), model.grid.Nx, 1)
    zC = repeat(reshape(model.grid.zC, 1, model.grid.Nz), model.grid.Nx, 1)

    u, v, w = model.velocities
    T, S = model.tracers

    sca(axs[2]); cla()
    title("Kinetic energy")
    pcolormesh(xC, zC, kinetic_energy_slice(data(u), data(v)))
    xlabel("\$ x \$ (m)")

    sca(axs[2]); cla()
    title("Vertical velocity")
    pcolormesh(xC, zF, data(w)[:, jhalf, :])
    xlabel("\$ x \$ (m)"); ylabel("\$ z \$ (m)")

    sca(axs[3]); cla()
    title("Salinity")
    pcolormesh(xC, zC, data(S)[:, jhalf, :])
    xlabel("\$ x \$ (m)")

    [ax.set_aspect(1) for ax in axs]
    pause(0.01)

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

    model.architecture == CPU() && makeplot!(axs, model)
end
