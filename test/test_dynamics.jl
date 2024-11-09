# test/test_dynamics.jl

using Test
using YAML

include("../src/MainModule.jl")
using .MainModule

# Load parameters from the config file
config_path = joinpath(@__DIR__, "..", "configs", "config.yaml")
params = Parameters.load_parameters(config_path)

# Create an instance of the RocketDynamics model
dynamics_model = RocketDynamics3dof(params)

# Test the dynamics function
@testset "RocketDynamics Tests" begin
    println("Testing dynamics function...")
    
    # Select a test state and control input
    x_test = params["x0"]
    u_test = params["u_guess"]
    
    # Call the dynamics function
    dx = dynamics(dynamics_model, x_test, u_test, params)
    
    # Check that the output is of correct size and type
    @test isa(dx, Vector{Float64}) || isa(dx, Vector{Float32})
    @test length(dx) == params["n_states"]
    
    # Print the output
    println("dx = ", dx)
    
    # Test the state Jacobian
    println("Testing state Jacobian...")
    A = state_jacobian(dynamics_model, x_test, u_test, params)
    @test size(A) == (params["n_states"], params["n_states"])
    
    # Print the state Jacobian
    println("A = ", A)
    
    # Test the control Jacobian
    println("Testing control Jacobian...")
    B = control_jacobian(dynamics_model, x_test, u_test, params)
    @test size(B) == (params["n_states"], params["n_controls"])
    
    # Print the control Jacobian
    println("B = ", B)
    
    # Test trajectory initialization
    println("Testing trajectory initialization...")
    x_init, u_init = initialize_trajectory(dynamics_model, params)
    @test size(x_init) == (params["N"], params["n_states"])
    @test size(u_init) == (params["N"] - 1, params["n_controls"])
    
    # Print the initial trajectories
    println("x_init[1] = ", x_init[1, :])
    println("u_init[1] = ", u_init[1, :])
end
