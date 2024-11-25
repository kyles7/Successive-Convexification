# test/test_convex_subproblem.jl
module TestConvexSubproblem

export runTestConvexSubproblem

# Change the working directory to the script's directory
cd(@__DIR__)


using ..MainModule
using LinearAlgebra
# using .Parameters
# using .AbstractDynamicsModel
# using .RocketDynamics
# using .Linearization
# using .ConvexSubproblemSolver
using Test


function runTestConvexSubproblem()
    # Load parameters
    config_path = joinpath(@__DIR__, "..", "..", "configs", "config6dof.yaml")
    params = Parameters.load_parameters(config_path)
    sigma = params["sigma"]
    # Initialize reference trajectories
    nondimensionalize!(params)
    X, U = initialize_trajectory6dof(params)
    X = Array{Float64,2}(X)
    U = Array{Float64,2}(U)

    # Calculate the discretization matrices
    A_list, B_list, C_list, S_list, Z_list = calculate_discretization(X, U, params["sigma"], params)

    # Solve the convex subproblem
    x_opt, u_opt = solve_convex_subproblem(A_list, B_list, C_list, S_list, Z_list, X, U, X, U, sigma, sigma, params)
    println("Convex subproblem solved successfully.")

    redimensionalize!(params)
    r_scale = norm(params["x0"][2:4])
    m_scale = params["m_wet"]
    
    # redimensionalize the trajectories
    redim_trajectory!(x_opt, u_opt, params)
    redim_trajectory!(X, U, params)
    
    # Print the results
    # println("Initial trajectory: ")
    # println("---------------------------------")
    # println("Reference state trajectory: ", X)
    # println("---------------------------------")
    # println("Reference control trajectory: ", U)
    # println("---------------------------------")
    # println("Results after redimensionalization: ")
    # println("---------------------------------")
    # println("Optimized state trajectory: ", x_opt)
    # println("---------------------------------")
    # println("Optimized control trajectory: ", u_opt)

    thrust_scale = 0.0002
    attitude_scale = 100

    println("Plotting...")
    # Currently, one of the plots can be displayed at a time, so uncomment only one

    # plot_trajectory(X, U, thrust_scale, attitude_scale) # Plot the initial trajectory
    # plot_trajectory(x_opt, u_opt, thrust_scale, attitude_scale) # Plot the optimized trajectory
    plot_thrust_magnitude(u_opt, params) # Plot the thrust magnitude over time


    # Test the results
    # @testset "Convex Subproblem Solver Tests" begin
    #     # Check dimensions of the optimized trajectories
    #     @test size(x_opt) == size(x_ref)
    #     @test size(u_opt) == size(u_ref)

    #     # Verify that the optimized trajectories satisfy the initial condition
    #     @test x_opt[1, :] ≈ params["x0"] atol=1e-6

    #     # Verify that the optimized trajectories satisfy the final conditioms
    #     @test x_opt[end, 1:3] ≈ params["rf"] atol=1e-6
    #     @test x_opt[end, 4:6] ≈ params["vf"] atol=1e-6
    #     @test x_opt[end, 7:10] ≈ params["qf"] atol=1e-6
    #     @test x_opt[end, 11:13] ≈ params["omegaf"] atol=1e-6
    #     # Optionally, you can test dynamics constraints by simulating the optimized trajectory
    #     # dt = params["dt"]
    #     # N = params["N"]
    #     # for k in 1:N-1
    #     #     xk = x_opt[k, :]
    #     #     uk = u_opt[k, :]
    #     #     xk1_expected = x_opt[k+1, :]

    #     #     # Compute the next state using the dynamics
    #     #     dx = dynamics6dof(xk, uk, params)
    #     #     xk1_simulated = xk + dt * dx

    #     #     # Compare the simulated next state with the optimized next state
    #     #     @test norm(xk1_expected - xk1_simulated) / norm(xk1_expected) ≤ 2e-3
    #     # end

    #     println("All tests passed for the convex subproblem solver.")
    # end
end # function runTestConvexSubproblem

end # module TestConvexSubproblem