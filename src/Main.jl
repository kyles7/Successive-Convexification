# src/Main.jl

module Main

using .DynamicsModels.RocketDynamics
using .Optimization.SuccessiveConvexification
using .Utilities.Parameters
using .Utilities.TrajectoryVisualization

function main()
    # Load parameters
    config_path = joinpath(@__DIR__, "..", "configs", "rocket_config.yaml")
    params = load_parameters(config_path)
    
    # Create an instance of RocketDynamics
    dynamics_model = RocketDynamics(params)
    
    # Initialize trajectories
    x_init, u_init = initialize_trajectory(dynamics_model, params)
    
    # Run the successive convexification algorithm
    x_opt, u_opt = successive_convexification(dynamics_model, x_init, u_init, params)
    
    # Visualize the optimized trajectory
    visualize_trajectory(x_opt, u_opt, params)
end

end # module Main
