# # Ocean convection example
#
# In this example, two-dimensional convection into a stratified fluid
# mixes a phytoplankton-like tracer. This example demonstrates how
#
#   * to set boundary conditions;
#   * to defined and insert a user-defined forcing function into a simulation.
#   * to use the `TimeStepWizard` to manage and adapt the simulation time-step.
# 
# To begin, we load Oceananigans, a plotting package, and a few miscellaneous useful packages.

using Oceananigans, PyPlot, Random, Printf

# ## Parameters
#
# We choose a modest two-dimensional resolution of 128² in a 64² m² domain ,
# implying a resolution of 0.5 m. Our fluid is initially stratified with 
# a squared buoyancy frequency
#
# $$ N\^2 = 10^{-5} \, \mathrm{s^{-2}} $$
#
# and a surface buoyancy flux 
#
# $$ Q_b2 = 10^{-8} \, \mathrm{m^3 \, s^{-2}} $$
#
# Because we use the physics-based convection whereby buoyancy flux by a
# positive vertical velocity implies positive flux, a positive buoyancy flux
# at the top of the domain carries buoyancy out of the fluid and causes convection.
# Finally, we end the simulation after 1 day.

Nz = 64
Lz = 64.0
N² = 1e-5
Qb = 1e-8
end_time = 1day
grid = RegularCartesianGrid(N = (Nz, 1, Nz), L = (Lz, Lz, Lz))
arch = CPU()

tracer_names = (:b, :nutrients, :plankton)
tracers = Oceananigans.TracerFields(arch, grid, tracer_names)
b, N, P = tracers

 uptake_half_saturation = k = 0.5
light_penetration_depth = δ = 16.0
         metabolic_loss = r = 0.07
   suraface_growth_rate = μ = 1 / hour

light_penetration(x, y, z) = μ * exp(z / δ)

nutrients_rhs  =  - (light_penetration * N / (k + N) - r) * P 
plankton_rhs   =    (light_penetration * N / (k + N) - r) * P

# ## Creating boundary conditions
#
# Create boundary conditions. Note that temperature is buoyancy in our problem.
#

buoyancy_bcs = HorizontallyPeriodicBCs(   top = BoundaryCondition(Flux, Qb), 
                                       bottom = BoundaryCondition(Gradient, N²))

## Instantiate the model
model = Model(architecture = arch,
                      grid = grid,
                   closure = ConstantIsotropicDiffusivity(ν=1e-4, κ=1e-4),
                  coriolis = FPlane(f=1e-4), 
                   tracers = tracers,
                  buoyancy = BuoyancyTracer(),
                   forcing = ModelForcing(nutrients = OperationForcing(nutrients_rhs), 
                                          plankton = OperationForcing(plankton_rhs)), #, herbavores=herbavores_rhs),
       boundary_conditions = BoundaryConditions(b=buoyancy_bcs)
)

## Set initial condition. Initial velocity and salinity fluctuations needed for AMD.
Ξ(z) = randn() * z / Lz * (1 + z / Lz) # noise
b₀(x, y, z) = N² * z + N² * Lz * 1e-6 * Ξ(z)
set!(model, b=b₀, nutrients=1, plankton=1)

## A wizard for managing the simulation time-step.
wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=90.0)

## Create a plot
fig, axs = subplots(ncols=3, figsize=(14, 6))

xC = repeat(model.grid.xC, 1, model.grid.Nz)
zF = repeat(reshape(model.grid.zF[1:end-1], 1, model.grid.Nz), model.grid.Nx, 1)
zC = repeat(reshape(model.grid.zC, 1, model.grid.Nz), model.grid.Nx, 1)

## Run the model
while model.clock.time < end_time
    update_Δt!(wizard, model)
    walltime = @elapsed time_step!(model, 10, wizard.Δt)

    sca(axs[1]); cla()
    pcolormesh(xC, zF, model.velocities.w[:, 1, :])
    title("Vertical velocity")
    xlabel("\$ x \$ (m)")
    ylabel("\$ z \$ (m)")

    sca(axs[2]); cla()
    pcolormesh(xC, zC, model.tracers.plankton[:, 1, :])
    title("Phytoplankton concentration")
    xlabel("\$ x \$ (m)")
    axs[2].tick_params(left=false, labelleft=false, right=false, labelright=false)

    sca(axs[3]); cla()
    pcolormesh(xC, zC, model.tracers.nutrients[:, 1, :])
    title("Nutrients")
    xlabel("\$ x \$ (m)")
    axs[3].tick_params(left=false, labelleft=false)

    suptitle(@sprintf("\$ t = %.2f\$ hours", model.clock.time / hour))
    [ax.set_aspect(1) for ax in axs]
    gcf(); pause(0.01)
end
