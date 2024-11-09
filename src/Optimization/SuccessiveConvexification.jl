# SuccessiveConvexification.jl

module SuccessiveConvexification

using JuMP, Gurobi
using ..DynamicsModels
using ..Utilities

function successive_convexification(model::DynamicsModel, x_init, u_init, params)
    max_iters = params["max_iters"]
    tolerance = params["tolerance"]
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