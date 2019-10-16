# # Two dimensional turbulence example
#
# In this example, we initialize a random velocity field and observe its viscous,
# turbulent decay in a two-dimensional domain. This example demonstrates:
#
#   * How to run a model with no buoyancy equation or tracers;
#   * How to create user-defined fields
#   * How to use differentiation functions
#
# For this example, we need `PyPlot` for plotting and `Statistics` for setting up
# a random initial condition with zero mean velocity.

using Oceananigans, PyPlot, Statistics, Oceananigans.AbstractOperations

# In addition to importing plotting and statistics packages, we import
# some types and functions from `Oceananigans` that will aid in the calculation
# and visualization of voriticty.

using Oceananigans: Face, Cell
using Oceananigans.TurbulenceClosures: ∂x_faa, ∂y_afa

# We instantiate the model with a simple isotropic diffusivity

model = Model(
        grid = RegularCartesianGrid(N=(128, 128, 1), L=(2π, 2π, 2π)),
    buoyancy = nothing, 
     tracers = nothing,
     closure = ConstantIsotropicDiffusivity(ν=1e-3, κ=1e-3)
)

# Our initial condition randomizes `u` and `v`. We also ensure that both have
# zero mean for purely aesthetic reasons.

u₀ = rand(size(model.grid)...)
u₀ .-= mean(u₀) 

set!(model, u=u₀, v=u₀)

# Finally, we create the vorticity field for storing `u` and `v`, initialize a
# figure, and run the model forward
u, v, w = model.velocities
vorticity_operation = ∂x(v) - ∂y(u)

ω = Field(Face, Face, Cell, model.architecture, model.grid) 
vorticity_computation = Computation(vorticity_operation, ω)

close("all")
fig, ax = subplots()

for i = 1:10
    time_step!(model, Nt=100, Δt=1e-1)

    compute!(vorticity_computation)

    cla()
    imshow(ω[:, :, 1])
    ax.axis("off")
    pause(0.1)
end
