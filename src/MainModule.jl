# src/MainModule.jl
module MainModule
    ######## use this block for test_dynamics ########
    # include("Utilities/Parameters.jl")  
    # include("DynamicsModels/RocketDynamics3dof.jl")  # use this for test_dynamics
    # using .Parameters, .RocketDynamics3dof
    # export Parameters, RocketDynamics_3dof, dynamics3dof, state_jacobian3dof, control_jacobian3dof, initialize_trajectory3dof

    ######## use this block for test_6dofdynamics #######
    include("Utilities/Parameters.jl")
    include("Utilities/TrajectoryVisualization.jl")
    include("Utilities/VarPlotter.jl")
    include("DynamicsModels/RocketDynamics3dof.jl")  # use this for test_dynamics
    include("DynamicsModels/RocketDynamics6dof.jl")  # use this for test_6dofdynamics
    include("Optimization/ConvexSubproblemSolver2.jl")
    include("Optimization/SuccessiveConvexification.jl")
    using .Parameters, .TrajectoryVisualization, .VarPlotter, .RocketDynamics6dof, .RocketDynamics3dof, .ConvexSubproblemSolver, .SuccessiveConvexification
    export Parameters, plot_trajectory, plot_thrust_magnitude, plot_angular_velocity_magnitude, plot_vector_component, RocketDynamics_6dof, dynamics6dof, state_jacobian6dof, control_jacobian6dof, initialize_trajectory6dof, quaternion_to_rotation_matrix, RocketDynamics_3dof, dynamics3dof, state_jacobian3dof, control_jacobian3dof, initialize_trajectory3dof, solve_convex_subproblem, calculate_discretization, x_nondim!, u_nondim!, x_redim!, u_redim!, nondimensionalize!, redimensionalize!, redim_trajectory!, integrate_nonlinear_piecewise, get_nonlinear_cost, successive_convexification

    ######## use this block for test_convex_subproblems #######
#     include("Utilities/Parameters.jl")
#     include("DynamicsModels/RocketDynamics6dof.jl")  # use this for test_6dofdynamics
#     include("Optimization/ConvexSubproblemSolver.jl")
#     # include("DynamicsModels/AbstractDynamicsModel.jl")
#     using .Parameters, .RocketDynamics, .ConvexSubproblemSolver
#     export Parameters, RocketDynamics_6dof, ConvexSubproblemSolver, solve_convex_subproblem, initialize_trajectory, state_jacobian, control_jacobian
# # #   

# include("Utilities/Parameters.jl")  
#     # include("DynamicsModels/AbstractDynamicsModel.jl")
#     #include("DynamicsModels/RocketDynamics.jl")  # use this for test_dynamics
#     include("DynamicsModels/RocketDynamics6dof.jl")  # use this for test_6dofdynamics

#     using .Parameters, .RocketDynamics
#     # using .AbstractDynamicsModel
    
#     # Export everything from submodules
#     export Parameters, RocketDynamics3dof, dynamics, state_jacobian, control_jacobian, initialize_trajectory
end
