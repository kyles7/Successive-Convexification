# test/test_convex_subproblem.jl
module TestSCVX

export runTestSCVX

# Change the working directory to the script's directory
cd(@__DIR__)


using ..MainModule
using LinearAlgebra

using Test


function runTestSCVX()
    # Load parameters
    config_path = joinpath(@__DIR__, "..", "..", "configs", "config6dof.yaml")
    params = Parameters.load_parameters(config_path)
    sigma = params["sigma"]
    nondimensionalize!(params)
    # Initialize reference trajectories
    all_X, all_U, all_sigma = successive_convexification()

    # get last trajectory
    X = all_X[end]
    U = all_U[end]
    sigma = all_sigma[end]

    redimensionalize!(params)
    r_scale = norm(params["x0"][2:4])
    m_scale = params["m_wet"]
    # redimensionalize the trajectories
    redim_trajectory!(X, U, params)

    thrust_scale = 0.00004
    attitude_scale = 30

    println("Plotting...")
    # Currently, one of the plots can be displayed at a time, so uncomment only one

    # plot_trajectory(X, U, thrust_scale, attitude_scale) # Plot the initial trajectory
    plot_trajectory(X, U, thrust_scale, attitude_scale) # Plot the optimized trajectory
    #plot_thrust_magnitude(u_opt, params) # Plot the thrust magnitude over time
end # function runTestConvexSubproblem

end # module TestConvexSubproblem