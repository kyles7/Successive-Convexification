# src/MainModule.jl
module MainModule
    include("Utilities/Parameters.jl")
    include("Utilities/TrajectoryVisualization.jl")
    include("Utilities/VarPlotter.jl")
    include("Utilities/SystemWriter.jl")
    include("DynamicsModels/RocketDynamics6dof.jl")  # use this for test_6dofdynamics
    include("Optimization/ConvexSubproblemSolver.jl")
    include("Optimization/SuccessiveConvexification.jl")
    using .Parameters, .TrajectoryVisualization, .VarPlotter, .SystemWriter, .RocketDynamics6dof, .ConvexSubproblemSolver, .SuccessiveConvexification
    export Parameters, plot_trajectory, plot_thrust_magnitude, plot_angular_velocity_magnitude, plot_vector_component, write_trajectory_data, RocketDynamics_6dof, dynamics6dof, state_jacobian6dof, control_jacobian6dof, initialize_trajectory6dof, quaternion_to_rotation_matrix, solve_convex_subproblem, calculate_discretization, x_nondim!, u_nondim!, x_redim!, u_redim!, nondimensionalize!, redimensionalize!, redim_trajectory!, integrate_nonlinear_piecewise, get_nonlinear_cost, successive_convexification
end
