# # Simple diffusion example
#
# This script provides our simplest example of Oceananigans.jl functionality:
# the diffusion of a one-dimensional Gaussian. This example demonstrates
#
#   * how to load `Oceananigans.jl`;
#   * how to instantiate an `Oceananigans.jl` `Model`;
#   * how to set an initial condition with a function;
#   * how to time-step a model forward, and finally
#   * how to look at results.
#
# ## Using `Oceananigans.jl`
#
# To use `Oceananigans.jl` after it has been installed, we bring
# `Oceananigans.jl` functions and names into our 'namespace' by writing

using Oceananigans, Oceananigans.AbstractOperations, PyPlot, Printf
using Oceananigans: Face, Cell

# ## Instantiating and configuring a `Model`
#
# To begin using Oceananigans, we instantiate a `Model` by calling the 
# `Model` constructor:
grid = RegularCartesianGrid((1, 1, 128), (1, 1, 1))

tracer_names = (:b,)
tracers = Oceananigans.TracerFields(CPU(), grid, tracer_names)
b = tracers.b

weird_diffusivity = @at (Cell, Cell, Cell) ∂z( (tanh(∂z(b) / 0.1) + 1) / 2 * ∂z(b))

buoyancy_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, 0.1))

model = Model(
    grid = grid,
    closure = ConstantIsotropicDiffusivity(κ = 0.1),
    tracers = tracers,
    buoyancy = nothing,
    forcing = ModelForcing(b=OperationForcing(weird_diffusivity)),
    boundary_conditions = BoundaryConditions(b=buoyancy_bcs)
)

# The keyword arguments `grid` and `closure` indicate that 
# our model grid is Cartesian with uniform grid spacing, that our diffusive
# stress and tracer fluxes are determined by diffusion with a constant 
# diffusivity `κ` (note that we do not use viscosity in this example).

# Note that by default, a `Model` has no-flux boundary condition on all
# variables. Next, we set an initial condition on our "passive tracer", 
# temperature. Our objective is to observe the diffusion of a Gaussian.

## Build a Gaussian initial condition function with width `δ`:
δ = 0.1
bᵢ(x, y, z) = exp( -(z + 0.5)^2 / (2δ^2) )

## Set `model.tracers.T` to the function `Tᵢ`:
set!(model, b=bᵢ)

# ## Running your first `Model`
#
# Finally, we time-step the model forward using the function
# `time_step!`, with a time-step size that ensures numerical stability.

## Time-scale for diffusion across a grid cell
cell_diffusion_time_scale = model.grid.Δz^2 / 1.0

## The function `time_step!` executes `Nt` time steps with step size `Δt`
## using a second-order Adams-Bashforth method
time_step!(model, Nt = 1000, Δt = 0.1 * cell_diffusion_time_scale)

# ## Visualizing the results
#
# We use `PyPlot.jl` to look at the results.

## A convenient function for generating a label with the Current model time
tracer_label(model) = @sprintf("\$ t=%.3f \$", model.clock.time)

## Create a figure with `PyPlot.jl`
close("all")
fig, ax = subplots()
title("Diffusion of a Gaussian")
xlabel("Tracer concentration")
ylabel(L"z")

## Plot initial condition
plot(bᵢ.(0, 0, model.grid.zC), model.grid.zC, "--", label=L"t=0")

## Plot current solution
plot(model.tracers.b[1, 1, :], model.grid.zC, label=tracer_label(model))
legend()

# Interesting! Running the model even longer makes even more interesting results.

for i = 1:3
    time_step!(model, Nt = 1000, Δt = 0.001 * cell_diffusion_time_scale)
    plot(model.tracers.b[1, 1, :], model.grid.zC, label=tracer_label(model))
end

legend()
