# test/test_convex_subproblem.jl

# Change the working directory to the script's directory
cd(@__DIR__)

# Include necessary modules
# include("../src/Utilities/Parameters.jl")
# include("../src/DynamicsModels/AbstractDynamicsModel.jl")
# include("../src/DynamicsModels/RocketDynamics.jl")
# include("../src/Utilities/Linearization.jl")
# include("../src/Optimization/ConvexSubproblemSolver.jl")
include("../src/MainModule.jl")
using .MainModule
using LinearAlgebra
# using .Parameters
# using .AbstractDynamicsModel
# using .RocketDynamics
# using .Linearization
# using .ConvexSubproblemSolver
using Test

# Load parameters
config_path = joinpath(@__DIR__, "..", "configs", "config6dof.yaml")
params = Parameters.load_parameters(config_path)

# Set additional parameters
params["Q"] = Diagonal(repeat([1.0], params["n_states"]))
params["R"] = Diagonal(repeat([0.1], params["n_controls"]))
params["QN"] = Diagonal(repeat([10.0], params["n_states"]))
params["max_iterations"] = 10
params["convergence_tolerance"] = 1e-3

# Define control limits (replace F_max and M_max with appropriate values)
F_max = 1e5  # Maximum thrust force (N)
M_max = 1e4  # Maximum moment (N·m)
params["u_min"] = [-F_max, -F_max, -F_max, -M_max, -M_max, -M_max]
params["u_max"] = [F_max, F_max, F_max, M_max, M_max, M_max]

# Create dynamics model instance
dynamics_model = RocketDynamics_6dof(params)
# print the type of dynamics_model
println(typeof(dynamics_model))
params["dynamics_model"] = dynamics_model

# Initialize reference trajectories
x_ref, u_ref = initialize_trajectory(dynamics_model, params)
x_ref = Array{Float64,2}(x_ref)
u_ref = Array{Float64,2}(u_ref)

# Select a time step to test
k_test = 1  # You can choose any time step from 1 to N-1

# Linearize dynamics at the selected time step
xk = x_ref[k_test, :]
uk = u_ref[k_test, :]

# Ensure quaternion is normalized
q_norm = sqrt(sum(xk[7:10].^2))
xk[7:10] /= q_norm

# Compute Jacobians using the linearization module
A_list = Matrix{Float64}[]
B_list = Matrix{Float64}[]
for k in 1:params["N"] - 1
    xk = x_ref[k, :]
    uk = u_ref[k, :]

    # Normalize quaternion
    q_norm = sqrt(sum(xk[7:10].^2))
    xk[7:10] /= q_norm

    # Compute Jacobians
    A = state_jacobian(dynamics_model, xk, uk, params)
    B = control_jacobian(dynamics_model, xk, uk, params)
    push!(A_list, A)
    push!(B_list, B)
end

# Solve the convex subproblem
x_opt, u_opt = solve_convex_subproblem(A_list, B_list, x_ref, u_ref, params)

# Test the results
@testset "Convex Subproblem Solver Tests" begin
    # Check dimensions of the optimized trajectories
    @test size(x_opt) == size(x_ref)
    @test size(u_opt) == size(u_ref)

    # Verify that the optimized trajectories satisfy the initial condition
    @test x_opt[1, :] ≈ params["x0"] atol=1e-6

    # Optionally, you can test dynamics constraints by simulating the optimized trajectory
    dt = params["dt"]
    N = params["N"]
    for k in 1:N-1
        xk = x_opt[k, :]
        uk = u_opt[k, :]
        xk1_expected = x_opt[k+1, :]

        # Compute the next state using the dynamics
        dx = dynamics(dynamics_model, xk, uk, params)
        xk1_simulated = xk + dt * dx

        # Compare the simulated next state with the optimized next state
        @test xk1_expected ≈ xk1_simulated atol=1e-3
    end

    println("All tests passed for the convex subproblem solver.")
end
