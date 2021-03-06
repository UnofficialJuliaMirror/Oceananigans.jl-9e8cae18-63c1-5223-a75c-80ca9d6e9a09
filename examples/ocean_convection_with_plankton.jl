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
# $ N² = 10⁻⁵ \rm{s⁻²} $
#
# and a surface buoyancy flux
#
# $ Q_b = 10⁻⁸ \rm{m³ s⁻²} $
#
# Because we use the physics-based convection whereby buoyancy flux by a
# positive vertical velocity implies positive flux, a positive buoyancy flux
# at the top of the domain carries buoyancy out of the fluid and causes convection.
# Finally, we end the simulation after 1 day.

Nz = 128
Lz = 64.0
N² = 1e-5
Qb = 1e-8
end_time = day / 2

# ## Creating boundary conditions
#
# Create boundary conditions. Note that temperature is buoyancy in our problem.
#

buoyancy_bcs = HorizontallyPeriodicBCs(   top = BoundaryCondition(Flux, Qb),
                                       bottom = BoundaryCondition(Gradient, N²))

# ## Define a forcing function
#
# Our forcing function roughly corresponds to the growth of phytoplankton in light
# (with a penetration depth of 16 meters here), and death due to natural mortality
# at a rate of 1 phytoplankton unit per second.

growth_and_decay = SimpleForcing((x, y, z, t) -> exp(z/16) - 1)

## Instantiate the model
model = Model(
                   grid = RegularCartesianGrid(size = (Nz, 1, Nz), length = (Lz, Lz, Lz)),
                closure = ConstantIsotropicDiffusivity(ν=1e-4, κ=1e-4),
               coriolis = FPlane(f=1e-4),
                tracers = (:b, :plankton),
               buoyancy = BuoyancyTracer(),
                forcing = ModelForcing(plankton=growth_and_decay),
    boundary_conditions = BoundaryConditions(b=buoyancy_bcs)
)

# Set makeplot = true to live-update a plot of vertical velocity, temperature, and salinity
# as the simulation runs.

makeplot = false

## Set initial condition. Initial velocity and salinity fluctuations needed for AMD.
Ξ(z) = randn() * z / Lz * (1 + z / Lz) # noise
b₀(x, y, z) = N² * z + N² * Lz * 1e-6 * Ξ(z)
set!(model, b=b₀)

## A wizard for managing the simulation time-step.
wizard = TimeStepWizard(cfl=0.1, Δt=1.0, max_change=1.1, max_Δt=90.0)

## Create a plot
fig, axs = subplots(ncols=2, figsize=(10, 6))

"""
    makeplot!(axs, model)

Make side-by-side x-z slices of vertical velocity and plankton associated with `model` in `axs`.
"""
function makeplot!(axs, model)
    xC = repeat(model.grid.xC, 1, model.grid.Nz)
    zF = repeat(reshape(model.grid.zF[1:end-1], 1, model.grid.Nz), model.grid.Nx, 1)
    zC = repeat(reshape(model.grid.zC, 1, model.grid.Nz), model.grid.Nx, 1)

    sca(axs[1]); cla()
    ## Calling the Array() constructor here allows plots to be made for GPU models:
    pcolormesh(xC, zF, Array(interior(model.velocities.w))[:, 1, :])
    title("Vertical velocity")
    xlabel("\$ x \$ (m)")
    ylabel("\$ z \$ (m)")

    sca(axs[2]); cla()
    ## Calling the Array() constructor here allows plots to be made for GPU models:
    pcolormesh(xC, zC, Array(interior(model.tracers.plankton))[:, 1, :])
    title("Phytoplankton concentration")
    xlabel("\$ x \$ (m)")
    axs[2].tick_params(left=false, labelleft=false)

    suptitle(@sprintf("\$ t = %.2f\$ hours", model.clock.time / hour))
    [ax.set_aspect(1) for ax in axs]
    pause(0.01)
    gcf()
    return nothing
end

# Run the model:

while model.clock.time < end_time
    update_Δt!(wizard, model)
    walltime = @elapsed time_step!(model, 100, wizard.Δt)

    ## Print a progress message
    @printf("progress: %.1f %%, i: %04d, t: %s, Δt: %s, wall time: %s\n", 
            model.clock.time / end_time * 100, model.clock.iteration, 
            prettytime(model.clock.time), prettytime(wizard.Δt), prettytime(walltime))

    makeplot && makeplot!(axs, model)
end

# Plot the result at the end
makeplot!(axs, model)
gcf()
