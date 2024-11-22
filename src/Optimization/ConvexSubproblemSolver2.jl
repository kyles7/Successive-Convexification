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
function solve_convex_subproblem(A_bar, B_bar, C_bar, S_bar, Z_bar, X, U, X_last, U_last, sigma, sigma_last, params::Dict)
    N = params["N"]
    n_states = params["n_states"]
    n_controls = params["n_controls"]
    T_max = params["T_max"]
    m_dry = params["m_dry"]
    T_min = params["T_min"]
    omega_max = params["omega_max"]

    
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    # define decision variables: state, control, nu, sigma 
    @variable(model, x[1:n_states, 1:N])
    @variable(model, u[1:n_controls, 1:N])
    @variable(model, nu[1:n_states, 1:N-1])
    @variable(model, sigma >=0)
    @variable(model, s_prime[N,1]>=0)

    # ------------------- Constraints --------------------------- #
    # ------------------- BOUNDARY CONDITIONS ------------------- #
    @constraint(model, x[:, 1] .== params["x0"])
    # final position constraint
    @constraint(model, x[1:3, N] .== params["xf"][1:3])
    # final velocity constraint
    @constraint(model, x[4:6, N] .== params["xf"][4:6])
    # final quaternion constraint
    @constraint(model, x[7:10, N] .== params["xf"][7:10])
    # final angular velocity constraint
    @constraint(model, x[11:13, N] .== params["xf"][11:13])
    # no final mass constraint 
    # ------------------- STATE CONSTRAINTS ---------------------- #
    # TODO: mass must be higher than dry mass
    @constraint(model, x[14, :] .>= params["m_dry"])
    # TODO: SOCC Glideslope Constraint
    # TODO: SOCC Angle Constraint
    # TODO: Angular velocity constraint
    for k in 1:N
        @constraint(model, [omega_max; x[11:13,k]] in SecondOrderCone())
    end
    # ------------------- CONTROL CONSTRAINTS -------------------- #
    # TODO: Thrust magnitude upper bound
    for k in 1:N
        @constraint(model, [T_max; u[:,k]] in SecondOrderCone())
    end
    # TODO: Linearize thrust lower bound
    #TODO: gimbal angle constraint
    # ------------------- DYNAMICS CONSTRAINTS ------------------- #
    for k in 1:N-1
        @constraint(model, x[:, k+1] .== 
        reshape(A_bar[:, k], n_states, n_states) * x[:, k] 
        + reshape(B_bar[:, k], n_states, n_controls) * u[:, k] 
        + reshape(C_bar[:, k], n_states, n_controls) * u[:, k+1]
        + S_bar[:, k] * sigma
        + Z_bar[:, k]
        + nu[:, k])
    end
    # ------------------- TRUST REGION -------------------------- #
    du = u - U_last
    dx = x - X_last
    ds = sigma - sigma_last 
    # TODO: sum of norms less than max trust region radius
   # @constraint(model, norm(du) + norm(dx) + norm(ds) <= params["tr_radius"])

    # ------------------- OBJECTIVE FUNCTION --------------------- #
    @objective(model, Min, 
        params["weight_sigma"] * sigma +
        #params["weight_nu"] * norm(nu) +
        #params["weight_nu"] * norm(nu) + 
        1e5 * sum(s_prime)
    )

    optimize!(model)

    # Check for solver convergence
    status = termination_status(model)
    if status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED
        error("Optimization failed with status $status")
    end

    # Extract optimized trajectories
    x_opt = Array{Float64}(undef, n_states, N)
    u_opt = Array{Float64}(undef, n_controls, N)
    for k in 1:N
        x_opt[:, k] = value.(x[:, k])
    end
    for k in 1:N
        u_opt[:, k] = value.(u[:, k])
    end

    return x_opt, u_opt
end



end # module ConvexSubproblemSolver
