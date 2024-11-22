# SuccessiveConvexification.jl

module SuccessiveConvexification

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
    X = Array{Float64,2}(x_ref)
    U = Array{Float64,2}(u_ref)
    sigma = params["sigma"]
    max_iters = params["max_iters"]
    tolerance = params["tolerance"]

    #----------------------- Main Loop ---------------------#
    last_nonlinear_cost = none 
    converged = false 
    for it in 1:max_iters
        # get current time 
        t0_it = time()
        println("------- Iteration $it -------")
        t0_tm = time()
        # Let discritized matrices
        A_bar, B_bar, C_bar, S_bar, z_bar = calculate_discretization(X, U, sigma, params)
        println("Time to calculate discretization: ", time() - t0_tm)
        z
        X_last = X
        U_last = U
        sigma_last = sigma
        weight_nu = params["weight_nu"]
        weight_sigma = params["weight_sigma"]
        tr_radius = params["tr_radius"]

        while true
            # Solve the convex subproblem
            # X_new, U_new, sigma_new, error = solve_convex_subproblem(A_bar, B_bar, C_bar, S_bar, z_bar, X, U, X_last, U_last sigma, params)
            
            # # Compute the linearized cost
            # cost_linearized = compute_linearized_cost(X, U, X_new, U_new, sigma, A_bar, B_bar, C_bar, S_bar, z_bar, params)
            
            # # Compute the actual reduction
            # actual_reduction = last_nonlinear_cost - cost_new
            
            # # Compute the predicted reduction
            # predicted_reduction = last_nonlinear_cost - cost_linearized
            
            # # Compute the ratio
            # ratio = actual_reduction / predicted_reduction
            
            # Update the trust region radius
            # if ratio < 0.25
            #     tr_radius *= 0.5
            # elseif ratio > 0.75 && norm(U_new - U) < 0.9 * tr_radius
            #     tr_radius = min(2 * tr_radius, params["tr_max"])
            # end
            if last_nonlinear_cost == none
                last_nonlinear_cost = nonlinear_cost 
                X = X_new
                U = U_new
                sigma = sigma_new
                break
            end
            actual_change = last_nonlinear_cost - nonlinear_cost  # delta_J
            predicted_change = last_nonlinear_cost - linear_cost  # delta_L
            
            print("")
            print(format_line("Virtual Control Cost", linear_cost_dynamics))
            print(format_line("Constraint Cost", linear_cost_constraints))
            print("")
            print(format_line("Actual change", actual_change))
            print(format_line("Predicted change", predicted_change))
            print("")
            print(format_line("Final time", sigma))
            print("")
            # Check the ratio and update the trajectory
            if ratio > 0.1
                X = X_new
                U = U_new
                sigma = sigma_new
                last_nonlinear_cost = cost_new
                break
            end
        end
        # Check convergence
        if has_converged(X, X_new, tolerance)
            println("Converged in $it iterations.")
            converged = true
            break
        end
        
        # Update current trajectory
        X = X_new
        U = U_new
    end


    x_current = x_init
    u_current = u_init
    
    for iter in 1:max_iters
        # Linearize dynamics around current trajectory
        A_list, B_list, c_list = linearize_dynamics(model, x_current, u_current, params)
        
        # Set up and solve the convex optimization problem
        x_next, u_next = solve_convex_subproblem(A_list, B_list, c_list, x_current, u_current, params)
        
        # Check convergence
        if has_converged(x_current, x_next, tolerance)
            println("Converged in $iter iterations.")
            break
        end
        
        # Update current trajectory
        x_current = x_next
        u_current = u_next
    end
    
    return x_current, u_current
end

end # module SuccessiveConvexification