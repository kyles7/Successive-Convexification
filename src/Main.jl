module Main

using .DynamicsModels
using .Optimization
using .Utilities

function main()
    # Load parameters
    params = load_parameters("config.yaml")
    
    # Select dynamics model
    dynamics_model = RocketDynamics(params)
    # Alternatively, for a different model:
    # dynamics_model = QuadrotorDynamics(params)
    
    # Initialize trajectory
    #x_init, u_init = initialize_trajectory(dynamics_model, params)
    
    # Run successive convexification
   # x_opt, u_opt = successive_convexification(dynamics_model, x_init, u_init, params)
    
    # Visualize results
   # visualize_trajectory(x_opt, u_opt, params)
end

end # module Main