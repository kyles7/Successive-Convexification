# SuccessiveConvexification.jl

module SuccessiveConvexification

export successive_convexification
using Test
using YAML
using LinearAlgebra
using ..MainModule
using DifferentialEquations
using Plots


function successive_convexification()
    config_path = joinpath(@__DIR__, "..", "..", "configs", "config6dof.yaml")
    params = Parameters.load_parameters(config_path)

    #----------------------- Initialization ---------------------#
    nondimensionalize!(params)
    X, U = initialize_trajectory6dof(params)
    X = Array{Float64,2}(X)
    U = Array{Float64,2}(U)
    sigma = params["sigma"]
    all_X = [copy(X)]
    all_U = [copy(U)]
    all_sigma = [sigma]
    max_iters = params["max_iters"]
    #----------------------- Main Loop ---------------------#
    last_nonlinear_cost = nothing  
    converged = false 
    for it in 1:max_iters
        # get current time 
        t0_it = time()
        println("------- Iteration $it -------")
        t0_tm = time()
        # Let discritized matrices
        A_bar, B_bar, C_bar, S_bar, z_bar = calculate_discretization(X, U, sigma, params)
        println("Time to calculate discretization: ", time() - t0_tm)
        X_last = X
        U_last = U
        sigma_last = sigma

        while true
            # Solve the convex subproblem
            X_new, U_new, sigma_new, nu_new, sprime_new = solve_convex_subproblem(A_bar, B_bar, C_bar, S_bar, z_bar, X, U, X_last, U_last, sigma_last, sigma_last, params)
            
            X_nl = integrate_nonlinear_piecewise(X_new, U_new, sigma_new, params)
            linear_cost_dynamics = norm(nu_new,1)
            nonlinear_cost_dynamics = norm(X_new - X_nl, 1)

            linear_cost_constraints = sprime_new
            nonlinear_cost_constraints = get_nonlinear_cost(X_new, U_new, params)

            linear_cost = linear_cost_dynamics + linear_cost_constraints
            nonlinear_cost = nonlinear_cost_dynamics + nonlinear_cost_constraints

            if isnothing(last_nonlinear_cost)
                last_nonlinear_cost = nonlinear_cost 
                X = X_new
                U = U_new
                sigma = sigma_new
                break
            end
            actual_change = last_nonlinear_cost - nonlinear_cost  # delta_J
            predicted_change = last_nonlinear_cost - linear_cost  # delta_L
            
            print("")
            println("Virtual Control Cost: ", linear_cost_dynamics)
            println("Constraint Cost: ", linear_cost_constraints)
            print("")
            println("Actual change: ", actual_change)
            println("Predicted change: ", predicted_change)
            print("")
            println("Final time: ", sigma)
            print("")
            # Check the ratio and update the trajectory
            if abs(predicted_change) < 1e-4
                converged = true
                break
            else
                rho = actual_change / predicted_change
                if rho < params["rho_0"]
                    params["tr_radius"] /= params["alpha"]
                    println("Shrinking trust region to ", params["tr_radius"])
                else
                    X = X_new
                    U = U_new
                    sigma = sigma_new
                    println("Accepting solution")
                    
                    if rho < params["rho_1"]
                        params["tr_radius"] /= params["alpha"]
                        println("Decreasing trust region to ", params["tr_radius"])
                    elseif rho > params["rho_2"]
                        params["tr_radius"] *= params["beta"]
                        println("Increasing trust region to ", params["tr_radius"])
                    end
                    last_nonlinear_cost = nonlinear_cost
                    break
                end
            end
            println('-' ^ 80)
        end
        println("Time for iteration $it: ", time() - t0_it)
        
        push!(all_X, copy(X))
        push!(all_U, copy(U))
        push!(all_sigma, sigma)
        # Check convergence
        if converged
            println("Converged in $it iterations.")
            converged = true
            break
        end

    end
    if !converged
        println("Did not converge in $max_iters iterations.")
    end

    #extract last trajectory
    X = all_X[end]
    U = all_U[end]
    sigma = all_sigma[end]
    #redimensionalize the trajectories
    # redimensionalize!(params)
    # redim_trajectory!(X, U, params)

    #plot the trajectory
    println("Plotting...")
    thrust_scale = 0.0002
    attitude_scale = 100
    plot_trajectory(X, U, thrust_scale, attitude_scale)
    return all_X, all_U, all_sigma
end 
end # module SuccessiveConvexification