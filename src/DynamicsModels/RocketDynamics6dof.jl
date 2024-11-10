module RocketDynamics

include("AbstractDynamicsModel.jl")
include("../Utilities/Parameters.jl")

export RocketDynamics, dynamics, initialize_trajectory

struct RocketDynamics_6dof <: AbstractDynamicsModel.DynamicsModel
    params::Dict
end

function dynamics(model::RocketDynamics_6dof, x::Vector{Float64}, u::Vector{Float64}, params::Dict) :: Vector{Float64}
    # Extract state variables
    # Position in inertial frame
    rx, ry, rz = x[1:3]
    # Velocity in inertial frame
    vx, vy, vz = x[4:6]
    # Quaternion (attitude)
    q0, q1, q2, q3 = x[7:10]
    # Angular velocity in body frame
    ωx, ωy, ωz = x[11:13]
    # Mass
    m = x[14]

    # Extract control inputs
    # Thrust force in body frame
    Fbx, Fby, Fbz = u[1:3]
    # Thrust moment in body frame
    Mbx, Mby, Mbz = u[4:6]

    # Parameters
    g_inertial = params["gravity_vector"]  # Gravitational acceleration vector [gx, gy, gz]
    g0 = params["standard_gravity"]        # Standard gravity (9.80665 m/s²)
    I_sp = params["I_sp"]                  # Specific impulse (s)
    I_body = params["inertia_matrix"]      # Inertia matrix in body frame (3x3)

    # Normalize quaternion to avoid drift
    q_norm = sqrt(q0^2 + q1^2 + q2^2 + q3^2)
    q0, q1, q2, q3 = q0 / q_norm, q1 / q_norm, q2 / q_norm, q3 / q_norm

    # Compute rotation matrix from body to inertial frame
    R_b_to_i = quaternion_to_rotation_matrix(q0, q1, q2, q3)

    # Translational dynamics
    dr = [vx, vy, vz]
    F_b = [Fbx, Fby, Fbz]
    dv = (1 / m) * (R_b_to_i * F_b) + g_inertial

    # Rotational dynamics
    q = [q0, q1, q2, q3]
    ω_b = [ωx, ωy, ωz]
    dq = 0.5 * quaternion_product(q, [0.0; ω_b])

    # Angular velocity derivative
    M_b = [Mbx, Mby, Mbz]
    dω = I_body \ (M_b - cross(ω_b, I_body * ω_b))

    # Mass dynamics
    F_b_magnitude = norm(F_b)
    dm = -F_b_magnitude / (g0 * I_sp)

    # Concatenate derivatives
    dx = vcat(dr, dv, dq, dω, dm)

    return dx
end

# Quaternion to rotation matrix
function quaternion_to_rotation_matrix(q0, q1, q2, q3) :: Matrix{Float64}
    R = [
        q0^2 + q1^2 - q2^2 - q3^2 2*(q1*q2 - q0*q3) 2*(q1*q3 + q0*q2);
        2*(q1*q2 + q0*q3) q0^2 - q1^2 + q2^2 - q3^2 2*(q2*q3 - q0*q1);
        2*(q1*q3 - q0*q2) 2*(q2*q3 + q0*q1) q0^2 - q1^2 - q2^2 + q3^2
    ]
    return R
end

# Quaternion product
function quaternion_product(q::Vector{Float64}, p::Vector{Float64}) :: Vector{Float64}
    q0, q1, q2, q3 = q
    p0, p1, p2, p3 = p
    s = q0 * p0 - q1 * p1 - q2 * p2 - q3 * p3
    x = q0 * p1 + q1 * p0 + q2 * p3 - q3 * p2
    y = q0 * p2 - q1 * p3 + q2 * p0 + q3 * p1
    z = q0 * p3 + q1 * p2 - q2 * p1 + q3 * p0
    return [s, x, y, z]
end

function state_jacobian(model::RocketDynamics_6dof, x::Vector{Float64}, u::Vector{Float64}, params::Dict) :: Matrix{Float64}
    A = ForwardDiff.jacobian(x -> dynamics(model, x, u, params), x)
    return A
end

function control_jacobian(model::RocketDynamics_6dof, x::Vector{Float64}, u::Vector{Float64}, params::Dict) :: Matrix{Float64}
    B = ForwardDiff.jacobian(u -> dynamics(model, x, u, params), u)
    return B
end

function initialize_trajectory(model::RocketDynamics_6dof, params::Dict) :: Tuple{Array{Float64,2}, Array{Float64,2}}
    N = params["N"]
    n_states = params["n_states"]
    n_controls = params["n_controls"]

    x_init = zeros(N, n_states)
    u_init = zeros(N - 1, n_controls)

    x_init[1, :] = params["x0"]

    # Simple initial guess for control inputs
    for k in 1:N-1
        xk = x_init[k, :]
        uk = params["u_guess"]
        dx = dynamics(model, xk, uk, params)
        x_init[k+1, :] = xk + params["dt"] * dx
        u_init[k, :] = uk
    end

    return x_init, u_init
end
    
end # module RocketDynamics6dof