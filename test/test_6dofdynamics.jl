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
params = load_parameters(config_path)

# Create dynamics model instance
dynamics_model = RocketDynamics_6dof(params)

# Test dynamics function
@testset "6DOF Rocket Dynamics Tests" begin
    x_test = params["x0"]
    u_test = params["u_guess"]

    dx = dynamics(dynamics_model, x_test, u_test, params)

    @test length(dx) == params["n_states"]

    println("dx = ", dx)
end
