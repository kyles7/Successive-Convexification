module ConvexSubproblemSolver

export solve_convex_subproblem

using JuMP
using Gurobi
using LinearAlgebra
using ..MainModule
using ForwardDiff
using ECOS
using Ipopt
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
function solve_convex_subproblem(A_bar, B_bar, C_bar, S_bar, Z_bar, X, U, X_last, U_last, sig, sigma_last, params::Dict)
    N = params["N"]
    n_states = params["n_states"]
    n_controls = params["n_controls"]
    T_max = params["T_max"]
    m_dry = params["m_dry"]
    T_min = params["T_min"]
    omega_max = params["rad_omega_max"]
    tan_gamma_gs = params["tan_gamma_gs"]
    cos_theta_max = params["cos_theta_max"]
    tan_delta_max = params["tan_delta_max"]

    
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    # define decision variables: state, control, nu, sigma 
    @variable(model, x[1:n_states, 1:N])
    @variable(model, u[1:n_controls, 1:N])
    @variable(model, nu[1:n_states, 1:N-1])
    @variable(model, sigma >=0)
    @variable(model, s_prime[1:N]>=0)

    # ------------------- Constraints --------------------------- #
    # ------------------- BOUNDARY CONDITIONS ------------------- #
    @constraint(model, x[:, 1] .== params["x0"])
    # final position constraint
    @constraint(model, x[2:end, N] .== params["xf"][2:end])
    # @constraint(model, x[1:3, N] .== params["xf"][1:3])
    # # final velocity constraint
    # @constraint(model, x[4:6, N] .== params["xf"][4:6])
    # # final quaternion constraint
    # @constraint(model, x[7:10, N] .== params["xf"][7:10])
    # # final angular velocity constraint
    # @constraint(model, x[11:13, N] .== params["xf"][11:13])
    # no final mass constraint 
    # ------------------- STATE CONSTRAINTS ---------------------- #
    # enforce non-negativity of z position
    #@constraint(model, x[4, :] .>= 0)
    # MIN  MASS CONSTRAINT
    @constraint(model, x[1, :] .>= params["m_dry"])

    # GLIDESLOPE CONSTRAINT
    # auxilarry variable t 
    @variable(model, t_k[1:N]>=0)
    for k in 1:N 
        @constraint(model, t_k[k] == x[4,k] / tan_gamma_gs)
        @constraint(model, [t_k[k]; x[2,k]; x[3,k]] in SecondOrderCone())
    end
    # TILT ANGLE CONSTRAINT
    c = sqrt((1 - cos_theta_max) / 2)
    for k in 1:N
        @constraint(model, [c; x[10,k]; x[11,k]] in SecondOrderCone())
    end
    # MAX ANGULAR VELOCITY CONSTRAINT
    for k in 1:N
        @constraint(model, [omega_max; x[12,k]; x[13,k]; x[14,k]] in SecondOrderCone())
    end
    # ------------------- CONTROL CONSTRAINTS -------------------- #
    # THRUST UPPER BOUND
    for k in 1:N
        @constraint(model, [T_max; u[:,k]] in SecondOrderCone())
    end
    # LOWER THRUST BOUND TODO: kinda sus, check this
    eps = 1e-6 
    u_last_p_unit = zeros(n_controls, N)
    for k in 1:N
        u_last_p_k = U_last[:, k]
        norm_u_last_p_k = norm(u_last_p_k) + eps
        u_last_p_unit[:, k] = u_last_p_k / norm_u_last_p_k
    end
    for k in 1:N
        lhs_k = sum(u_last_p_unit[:, k] .* u[:, k])  # scalar projecton
        @constraint(model, T_min - lhs_k <= s_prime[k])
    end
    # GIMBAL ANGLE CONSTRAINT
    # @variable(model, t_u[1:N]>=0) # aux var 
    # for k in 1:N
    #     @constraint(model, t_u[k] == tan_delta_max * u[3, k])
    #     @constraint(model, [t_u[k]; u[1,k]; u[2, k]] in SecondOrderCone())
    # end
    # # GIMBAL ANGLE CONSTRAINT
    @variable(model, t_u[1:N]>=0) # aux var 
    for k in 1:N
        @constraint(model, t_u[k] == 2*tan_delta_max * u[3, k])
        @constraint(model, [t_u[k]; u[1:2,k]] in SecondOrderCone())
    end

    # @variable(model, t_u[1:N]>=0) # aux var 
    # for k in 1:N
    #     #@constraint(model, t_u[k] == 2*tan_delta_max * u[3, k])
    #     @constraint(model, [2*tan_delta_max * u[3, k]; u[1,k]; u[2,k]] in SecondOrderCone())
    #    # @constraint(model, u[1,k]^2 + u[2,k]^2 <= 2 * tan_delta_max * u[3,k]^2)
    # end

    
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
    # TRUST REGION CONSTRAINT
    # @variable(model, dx_abs[1:n_states, 1:N] >=0)
    # @variable(model, du_abs[1:n_controls, 1:N] >=0)
    # @variable(model, ds_abs >=0)

    # for i in 1:n_states
    #     for k in 1:N
    #         @constraint(model, dx_abs[i,k] >= dx[i,k])
    #         @constraint(model, dx_abs[i,k] >= -dx[i,k])
    #     end
    # end
    # for i in 1:n_controls
    #     for k in 1:N
    #         @constraint(model, du_abs[i,k] >= du[i,k])
    #         @constraint(model, du_abs[i,k] >= -du[i,k])
    #     end
    # end
    # @constraint(model, ds_abs >= ds)
    # @constraint(model, ds_abs >= -ds)

    # norm_dx = sum(dx_abs)
    # norm_du = sum(du_abs)
    # norm_ds = ds_abs
    
    #Try 2: compute norms via dot products
    # for k in 1:N
    #     du_k_norm = du[:,k] * du[:,k]
    #     dx_k_norm = dx[:,k] * dx[:,k]
    #     @constraint(model, du_k_norm <= params["tr_radius"])
    #     @constraint(model, dx_k_norm <= params["tr_radius"])
    # end

    # Try 2: 2 norms
    @constraint(model, [params["tr_radius"]; vec(dx); vec(du); ds] in SecondOrderCone())
    # ds_norm = ds * ds
    # @constraint(model, ds_norm <= params["tr_radius"])
    # @constraint(model, norm_dx + norm_du + norm_ds <= params["tr_radius"])
#    @constraint(model, norm(du) + norm(dx) + norm(ds) <= params["tr_radius"])
#     @constraint(model, norm_ds^2 + norm_du^2 + norm_dx^2 <= params["tr_radius"]^2)
    # ------------------- OBJECTIVE FUNCTION --------------------- #
    @variable(model, nu_abs[1:n_states, 1:N-1] >=0)
    for i in 1:n_states
        for k in 1:N-1
            @constraint(model, nu_abs[i,k] >= nu[i,k])
            @constraint(model, nu_abs[i,k] >= -nu[i,k])
        end
    end
    norm_nu = sum(nu_abs)
    sum_s_prime = sum(s_prime)
    @objective(model, Min, 
        params["weight_sigma"] * sigma +
        params["weight_nu"] * norm_nu +
        1e5 * sum_s_prime
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
    sigma_new = value(sigma)
    nu_new = Array{Float64}(undef, n_states, N-1)
    for k in 1:N-1
        nu_new[:, k] = value.(nu[:, k])
    end

    sprime_new = sum(value.(s_prime))
    return x_opt, u_opt, sigma_new, nu_new, sprime_new  
end



end # module ConvexSubproblemSolver
