using Oceananigans, PyPlot, Random, Printf

using Oceananigans: Face, Cell

# Resolution
Nx = 128
Ny = 128
Nz = 64 

# Domain size
Lx = 1000e3
Ly = 1000e3
Lz = 800

# Physical parameters
 f = 1e-4     # Coriolis parameter
 α = 0.0      # Thermal expansion coefficient
 β = 8e-4     # Haline contraction coefficient

surface_salinity = 10
 bottom_salinity = 20

const τ₀ = 1e-5     # Wind stress magnitude
const Δτ = 70e3     # Width of wind stress

# Simulation end time
end_time = 30day

# Wind stress / velocity flux forcing functions
@inline τˣ(x, y, t) = τ₀ * y/Δτ * exp(-(x^2 + y^2) / 2Δτ^2)
@inline τʸ(x, y, t) = τ₀ * x/Δτ * exp(-(x^2 + y^2) / 2Δτ^2)

u_bcs = HorizontallyPeriodicBCs(   top = FunctionBoundaryCondition(Flux, :z, Face, Cell, τˣ),
                                bottom = BoundaryCondition(Value, 0))

v_bcs = HorizontallyPeriodicBCs(   top = FunctionBoundaryCondition(Flux, :z, Cell, Face, τʸ),
                                bottom = BoundaryCondition(Value, 0))

S_bcs = HorizontallyPeriodicBCs(   top = BoundaryCondition(Value, surface_salinity),
                                bottom = BoundaryCondition(Value, bottom_salinity))

model = Model(
         architecture = CPU(),
                 grid = RegularCartesianGrid(size=(Nx, Ny, Nz), x=(-Lx/2, Lx/2), y=(-Ly/2, Ly/2), z=(-Lz, 0)),
             coriolis = FPlane(f=f),
             buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=α, β=β)),
              closure = ConstantAnisotropicDiffusivity(νv=1e-2, νh=10, κv=1e-2, κh=10),
  boundary_conditions = BoundaryConditions(u=u_bcs, v=v_bcs, S=S_bcs)
)

# Initial condition
h = 100
initial_salinity(x, y, z) = (surface_salinity - bottom_salinity) * exp(-z^2 / 2h^2) + bottom_salinity

set!(model, S=initial_salinity)

wizard = TimeStepWizard(cfl=0.005, Δt=minute, max_change=1.1, max_Δt=20minute)

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
    title("\$ u(x, y, z=0) \$")
    imshow(interior(model.velocities.u)[:, :, Nz])

    kplot = Nz - 2
    sca(axs[2]); cla()
    title("\$ w(x, y, z=$(model.grid.zF[kplot])")
    imshow(interior(model.velocities.w)[:, :, kplot])
    
    axs[1].axes("off")
    axs[2].axes("off")
end
