using Oceananigans, PyPlot, Random, Printf

#using Oceananigans.AbstractOperations
#using Oceananigans: Face, Cell

Nx = 128
Ny = 128
Nz = 16 

Lx = 1000e3
Ly = 1000e3
Lz = 800

end_time = 30day

 f = 1e-4
 α = 0.0      # Thermal expansion coefficient
 β = 8e-4     # Haline contraction coefficient

surface_salinity = 10
 bottom_salinity = 20

const τ₀ = 1e-4     # Kinematic wind stress magnitude
const Δτ = 50e3     # Kinematic wind stress magnitude

@inline τˣ(x, y, t) = τ₀ * y/Δτ * exp(-(x^2 + y^2) / 2Δτ^2)
@inline τʸ(x, y, t) = τ₀ * x/Δτ * exp(-(x^2 + y^2) / 2Δτ^2)

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
                 grid = RegularCartesianGrid(size=(Nx, Ny, Nz), x=(-Lx/2, Lx/2), y=(-Ly/2, Ly/2), z=(-Lz, 0)),
             coriolis = FPlane(f=f),
             buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=α, β=β)),
              closure = ConstantAnisotropicDiffusivity(νv=1e-2, νh=1000, κv=1e-2, κh=1000),
  boundary_conditions = BoundaryConditions(u=u_bcs, v=v_bcs, S=S_bcs)
)

## Temperature initial condition: a stable density tradient with random noise superposed.
h = 100
initial_salinity(x, y, z) = surface_salinity * (exp(-z^2 / 2h^2) - exp(-Lz^2 / 2h^2)) + bottom_salinity

set!(model, S=initial_salinity)

wizard = TimeStepWizard(cfl=0.01, Δt=minute, max_change=1.1, max_Δt=4hour)

# A diagnostic that returns the maximum absolute value of `w` by calling
# `wmax(model)`:

wmax = FieldMaximum(abs, model.velocities.w)

# We also create a figure and define a plotting function for live plotting of results.

fig, axs = subplots(ncols=2, figsize=(12, 6))

xF_xy = repeat(model.grid.xF[2:end], 1, model.grid.Ny)
yC_xy = repeat(reshape(model.grid.yC, 1, model.grid.Ny), model.grid.Nx, 1)

xC_xz = repeat(model.grid.xC, 1, model.grid.Nz)
zF_xz = repeat(reshape(model.grid.zF[2:end], 1, model.grid.Nz), model.grid.Nx, 1)
zC_xz = repeat(reshape(model.grid.zC, 1, model.grid.Nz), model.grid.Nx, 1)

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

    sca(axs[1]); cla()
    pcolormesh(xF_xy, yC_xy, interior(model.velocities.u)[:, :, Nz])

    sca(axs[2]); cla()
    pcolormesh(xC_xz, zC_xz, interior(model.tracers.S)[:, round(Int, Ny/2), :])
end
