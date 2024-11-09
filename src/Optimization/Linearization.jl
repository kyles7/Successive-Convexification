# Linearization.jl

module Linearization

using ..DynamicsModels

function linearize_dynamics(model::DynamicsModel, x_current, u_current, params)
    N = length(x_current)
    A_list = []
    B_list = []
    c_list = []
    dt = params["dt"]
    
    for k in 1:N-1
        xk = x_current[k]
        uk = u_current[k]
        
        # Compute dynamics and Jacobians
        fk = dynamics(model, xk, uk, params)
        Ak = state_jacobian(model, xk, uk, params)
        Bk = control_jacobian(model, xk, uk, params)
        
        # Discretize dynamics (Euler integration)
        A_k = I + dt * Ak
        B_k = dt * Bk
        c_k = dt * (fk - Ak * xk - Bk * uk)
        
        push!(A_list, A_k)
        push!(B_list, B_k)
        push!(c_list, c_k)
    end
    
    return A_list, B_list, c_list
end

end # module Linearization
