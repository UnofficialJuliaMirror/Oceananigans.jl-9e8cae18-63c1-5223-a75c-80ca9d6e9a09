using Oceananigans, PyPlot, Random, Printf, Oceananigans.AbstractOperations

using Oceananigans: Face, Cell

Nx = 64 
Ny = 64 
Nz = 32 

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
const Δτ = 100e3     # Kinematic wind stress magnitude

#τˣ = BoundaryConditionFunction{:z, Face, Cell}((x, y, t) -> τ₀ * y/Δτ * exp(-y^2 / 2Δτ^2))
#τʸ = BoundaryConditionFunction{:z, Cell, Face}((x, y, t) -> τ₀ * x/Δτ * exp(-x^2 / 2Δτ^2))

@inline τˣ(x, y, t) = τ₀ * y/Δτ * exp(-y^2 / 2Δτ^2)
@inline τʸ(x, y, t) = τ₀ * x/Δτ * exp(-x^2 / 2Δτ^2)

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

wizard = TimeStepWizard(cfl=0.2, Δt=minute, max_change=1.1, max_Δt=1hour)

# A diagnostic that returns the maximum absolute value of `w` by calling
# `wmax(model)`:

wmax = FieldMaximum(abs, model.velocities.w)

# We also create a figure and define a plotting function for live plotting of results.

fig, axs = subplots()

function xyslice(a) 
    xC = repeat(model.grid.xC, 1, model.grid.Ny)
    yC = repeat(reshape(model.grid.yC, 1, model.grid.Ny), model.grid.Nx, 1)
    pcolormesh(xC, zC, data(a)[:, :, model.grid.Nz])
    return nothing
end

u, v, w = model.velocities
kinetic_energy = Computation(@at (Cell, Cell, Cell) u^2 + v^2, model)

"""
    makeplot!(axs, model)

Make a triptych of x-z slices of vertical velocity, temperature, and salinity
associated with `model` in `axs`.
"""
function makeplot!(axs, model)
    ke = kinetic_energy(model)
    xC = repeat(model.grid.xC, 1, model.grid.Ny)
    yC = repeat(reshape(model.grid.yC, 1, model.grid.Ny), model.grid.Nx, 1)

    pcolormesh(xC, zC, data(ke)[:, :, model.grid.Nz])

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
