module ConvexSubproblemSolver

export solve_convex_subproblem

using JuMP
using Gurobi
using LinearAlgebra
using ..MainModule
using ForwardDiff
using ECOS
using Ipopt

"""
    solve_convex_subproblem(A_list, B_list, x_ref, u_ref, params) -> (Array, Array)

Solves the convex subproblem in the Scvx algorithm using Gurobi.

# Arguments
- A_bar::Array: A matrix.
- B_bar::Array: B matrix.
- C_bar::Array: C matrix.
- S_bar::Array: S matrix.
- Z_bar::Array: Z matrix.
- X::Array: State trajectory.
- U::Array: Control trajectory.
- X_last::Array: Last state trajectory.
- U_last::Array: Last control trajectory.
- sig::Float64: Time guess.
- sigma_last::Float64: Last time Guess.

# Returns
- `(x_opt::Array, u_opt::Array)`: Optimized state and control trajectories.
- `sigma_new::Float64`: Optimized time.
- `nu_new::Array`: Optimized virtual control variables.
- `sprime_new::Float64`: Optimized slack variable.
"""
function solve_convex_subproblem(A_bar, B_bar, C_bar, S_bar, Z_bar, X_last, U_last, sigma_last, params::Dict)
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
    # ------------------- STATE CONSTRAINTS ---------------------- #
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
    @variable(model, t_u[1:N]>=0) # aux var 
    for k in 1:N
        @constraint(model, t_u[k] == 2*tan_delta_max * u[3, k])
        @constraint(model, [t_u[k]; u[1:2,k]] in SecondOrderCone())
    end
    
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
    @constraint(model, [params["tr_radius"]; vec(dx); vec(du); ds] in SecondOrderCone())
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
