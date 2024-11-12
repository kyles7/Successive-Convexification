# test/test_scvx.jl

include("../src/Utilities/Parameters.jl")
include("../src/DynamicsModels/AbstractDynamicsModel.jl")
include("../src/DynamicsModels/RocketDynamics.jl")
include("../src/Utilities/Linearization.jl")
include("../src/Optimization/ConvexSubproblemSolver.jl")

using .Parameters
using .AbstractDynamicsModel
using .RocketDynamics
using .Linearization
using .ConvexSubproblemSolver

# Load parameters
config_path = joinpath(@__DIR__, "..", "config.yaml")
params = load_parameters(config_path)

# Set additional parameters
params["Q"] = Diagonal(repeat([1.0], params["n_states"]))
params["R"] = Diagonal(repeat([0.1], params["n_controls"]))
params["QN"] = Diagonal(repeat([10.0], params["n_states"]))
params["max_iterations"] = 20
params["convergence_tolerance"] = 1e-3

# Define control limits
F_max = 1e5  # Replace with appropriate value
M_max = 1e4  # Replace with appropriate value
params["u_min"] = [-F_max, -F_max, -F_max, -M_max, -M_max, -M_max]
params["u_max"] = [F_max, F_max, F_max, M_max, M_max, M_max]

# Create dynamics model instance
dynamics_model = RocketDynamics(params)
params["dynamics_model"] = dynamics_model

# Run Scvx algorithm
x_opt, u_opt = scvx_algorithm(params)

# Print results
println("Optimized state trajectory:")
display(x_opt)
println("Optimized control trajectory:")
display(u_opt)
