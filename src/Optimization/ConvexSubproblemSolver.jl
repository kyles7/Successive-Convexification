module ConvexSubproblemSolver

export solve_convex_subproblem

using JuMP
using Gurobi
using LinearAlgebra
using ..MainModule
using ForwardDiff
#using ..AbstractDynamicsModel
# include("../DynamicsModels/AbstractDynamicsModel.jl")
# using .AbstractDynamicsModel

"""
    solve_convex_subproblem(A_list, B_list, x_ref, u_ref, params) -> (Array, Array)

Solves the convex subproblem in the Scvx algorithm using Gurobi.

# Arguments
- `A_list::Vector{Matrix}`: List of state Jacobian matrices at each time step.
- `B_list::Vector{Matrix}`: List of control Jacobian matrices at each time step.
- `x_ref::Array`: Reference state trajectory (N x n_states).
- `u_ref::Array`: Reference control trajectory (N-1 x n_controls).
- `params::Dict`: Dictionary containing problem parameters.

# Returns
- `(x_opt::Array, u_opt::Array)`: Optimized state and control trajectories.
"""
function solve_convex_subproblem(A_list::Vector{<:AbstractMatrix}, B_list::Vector{<:AbstractMatrix}, x_ref::AbstractArray, u_ref::AbstractArray, params::Dict)
    N = params["N"]
    n_states = params["n_states"]
    n_controls = params["n_controls"]
    dt = params["dt"]
    m_dry = params["m_dry"]

    # Extract weighting matrices for the objective function
    Q = params["Q"]        # State weighting matrix
    R = params["R"]        # Control weighting matrix
    QN = params["QN"]      # Terminal state weighting matrix (if any)

    # Create a JuMP model with Gurobi optimizer
    model = Model(Gurobi.Optimizer)

    # Suppress output from the solver (optional)
    set_silent(model)

    # Define decision variables
    @variable(model, x[1:N, 1:n_states])
    @variable(model, u[1:N-1, 1:n_controls])

    # Initial condition constraint
    @constraint(model, x[1, :] .== params["x0"])

    # Dynamics constraints
    for k in 1:N-1
        A = A_list[k]
        B = B_list[k]
        xk_ref = x_ref[k, :]
        uk_ref = u_ref[k, :]
        fk = dynamics_residual(xk_ref, uk_ref, params)

        @constraint(model, x[k+1, :] .== x[k, :] + dt * (A * (x[k, :] - xk_ref) + B * (u[k, :] - uk_ref) + fk))
    end

    # Objective function
    @expression(model, obj, 0)
    for k in 1:N-1
        @expression(model, obj += (x[k, :] - x_ref[k, :])' * Q * (x[k, :] - x_ref[k, :]) + (u[k, :] - u_ref[k, :])' * R * (u[k, :] - u_ref[k, :]))
    end
    # Terminal cost
    @expression(model, obj += (x[N, :] - x_ref[N, :])' * QN * (x[N, :] - x_ref[N, :]))

    @objective(model, Min, obj)

    # Control constraints (if any)
    if haskey(params, "u_min") && haskey(params, "u_max")
        u_min = params["u_min"]
        u_max = params["u_max"]
        @constraint(model, [k=1:N-1, j=1:n_controls], u_min[j] <= u[k, j] <= u_max[j])
    end

    # State constraints (if any)
    # mass constraint - mass must stay above dry mass
    mass_index = params["mass_index"]
    @constraint(model, [k=1:N], x[k, mass_index] >= 0)  # Altitude constraint
    #TODO: implement thrust constraint
    #TODO: implement tilt constraint

    # Solve the optimization problem
    optimize!(model)

    # Check for solver convergence
    status = termination_status(model)
    if status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED
        error("Optimization failed with status $status")
    end

    # Extract optimized trajectories
    x_opt = Array{Float64}(undef, N, n_states)
    u_opt = Array{Float64}(undef, N-1, n_controls)
    for k in 1:N
        x_opt[k, :] = value.(x[k, :])
    end
    for k in 1:N-1
        u_opt[k, :] = value.(u[k, :])
    end

    return x_opt, u_opt
end

# Helper function to compute the dynamics residual
function dynamics_residual(xk::AbstractVector, uk::AbstractVector, params::Dict) :: AbstractVector
    f = dynamics6dof(xk, uk, params)
    # TODO: replace forward diffs with jacobian functions
    # A = ForwardDiff.jacobian(x_state -> dynamics6dof(x_state, uk, params), xk)
    # B = ForwardDiff.jacobian(u_state -> dynamics6dof(xk, u_state, params), uk)
    A = state_jacobian6dof(xk, uk, params)
    B = control_jacobian6dof(xk, uk, params)
    residual = f - A * xk - B * uk
    return residual
end


end # module ConvexSubproblemSolver
