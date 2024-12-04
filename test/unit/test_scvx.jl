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

    # print the final sigma (time)
    println("Sigma = ", sigma)
    redimensionalize!(params)
    r_scale = norm(params["x0"][2:4])
    m_scale = params["m_wet"]
    # redimensionalize all the trajectories
    for i in 1:length(all_X)
        redim_trajectory!(all_X[i], all_U[i], params)
    end

    thrust_scale = 0.00004
    attitude_scale = 30

    # Plotting now disabled by default
    # println("Plotting...")
    # Currently, one of the plots can be displayed at a time, so uncomment only one
    # plot_trajectory(X, U, thrust_scale, attitude_scale) # Plot the initial trajectory
    # plot_trajectory(X, U, thrust_scale, attitude_scale) # Plot the optimized trajectory
    println("Writing trajectory data...")
    write_trajectory_data("trajectory_data.h5", all_X, all_U, all_sigma, params) # Write the trajectory data to a file
end # function runTestConvexSubproblem

end # module TestConvexSubproblem