# test/test_6dofdynamics.jl
using Test
using YAML
# include("../src/Utilities/Parameters.jl")
# include("../src/DynamicsModels/AbstractDynamicsModel.jl")
# include("../src/DynamicsModels/RocketDynamics6dof.jl")
include("../src/MainModule.jl")
using .MainModule


# Load parameters
config_path = joinpath(@__DIR__, "..", "configs", "config6dof.yaml")
params = Parameters.load_parameters(config_path)

# Create dynamics model instance
dynamics_model = RocketDynamics_6dof(params)

# Test dynamics function
@testset "6DOF Rocket Dynamics Tests" begin
    x_test = params["x0"]
    u_test = params["u_guess"]

    dx = dynamics(dynamics_model, x_test, u_test, params)

    @test length(dx) == params["n_states"]

    println("dx = ", dx)

    # Test quaternion to rotation matrix
    q0, q1, q2, q3 = x_test[7:10]
    R = quaternion_to_rotation_matrix(q0, q1, q2, q3)
    println("R = ", R)

    # Test the state Jacobian
    println("Testing state Jacobian...")
    A = state_jacobian(dynamics_model, x_test, u_test, params)
    @test size(A) == (params["n_states"], params["n_states"])
    # # Print the state Jacobian
    println("A = ", A)

    # Test the control Jacobian
    println("Testing control Jacobian...")
    B = control_jacobian(dynamics_model, x_test, u_test, params)
    @test size(B) == (params["n_states"], params["n_controls"])
    # # Print the control Jacobian
    println("B = ", B)
end
