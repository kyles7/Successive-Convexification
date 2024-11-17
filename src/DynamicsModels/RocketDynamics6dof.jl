module RocketDynamics6dof

using LinearAlgebra
using ForwardDiff

export RocketDynamics_6dof, dynamics6dof, state_jacobian6dof, control_jacobian6dof, initialize_trajectory6dof, quaternion_to_rotation_matrix

struct RocketDynamics_6dof 
    params::Dict
end

function dynamics6dof(x::AbstractVector, u::AbstractVector, params::Dict) :: AbstractVector{Real}
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
    # Position of engine wrt center of mass TODO: add to params
    rTB = [0.0; 0.0; -14.0]
    # Extract control inputs
    TBx, TBy, TBz = u[1:3]
    # Compute thrust and moment vectors
    FIx, FIy, FIz = quaternion_to_rotation_matrix(q0, q1, q2, q3) * [TBx; TBy; TBz]
    MBx, MBy, MBz = skew_symmetric3d(rTB) * [TBx; TBy; TBz]

    # # Thrust force in body frame
    # Fbx, Fby, Fbz = u[1:3]
    # # Thrust moment in body frame
    # Mbx, Mby, Mbz = u[4:6]

    #TODO: switch control vector to T

    # Parameters
    g_inertial = params["gravity_vector"]  # Gravitational acceleration vector [gx, gy, gz]
    g0 = params["standard_gravity"]        # Standard gravity (9.80665 m/s²)
    I_sp = params["I_sp"]                  # Specific impulse (s)
    I_body = hcat(params["inertia_matrix"]...)    # Inertia matrix in body frame (3x3)

    # Normalize quaternion to avoid drift
    q_norm = sqrt(q0^2 + q1^2 + q2^2 + q3^2)
    q0, q1, q2, q3 = q0 / q_norm, q1 / q_norm, q2 / q_norm, q3 / q_norm

    # Compute rotation matrix from body to inertial frame
    R_b_to_i = quaternion_to_rotation_matrix(q0, q1, q2, q3)

    # Translational dynamics
    dr = [vx, vy, vz]
    # F_b = [Fbx, Fby, Fbz]
    # dv = (1 / m) * (R_b_to_i * F_b) + g_inertial
    dv = (1 / m) * [FIx; FIy; FIz] + g_inertial 

    # Rotational dynamics
    q = [q0, q1, q2, q3]
    ω_b = [ωx, ωy, ωz]
    # dq = 0.5 * quaternion_product(q, [0.0; ω_b])
    dq = 0.5 * skew_symmetric4d(ω_b) * q
    # Angular velocity derivative
    M_b = [MBx, MBy, MBz]
    # dω = I_body \ (M_b - cross(ω_b, I_body * ω_b))
    dω = I_body \ (M_b - skew_symmetric3d(ω_b) * I_body * ω_b)
    # Mass dynamics
    # F_b_magnitude = norm(F_b)
    # dm = -F_b_magnitude / (g0 * I_sp)
    mag_TB = norm([TBx; TBy; TBz])
    dm = -mag_TB / (g0 * I_sp)
    # Concatenate derivatives
    dx = vcat(dr, dv, dq, dω, dm)

    return dx
end

# Quaternion to rotation matrix #TODO: move to utils, rename to DCM
function quaternion_to_rotation_matrix(q0::T, q1::T, q2::T, q3::T) :: AbstractMatrix{T} where T <: Real
    R = zeros(T, 3, 3)  # Generic matrix with type T
    R[1, 1] = 1 - 2 * (q2^2 + q3^2)
    R[1, 2] = 2 * (q1 * q2 - q3 * q0)
    R[1, 3] = 2 * (q1 * q3 - q2 * q0)
    R[2, 1] = 2 * (q1 * q2 - q3 * q0)
    R[2, 2] = 1 - 2 * (q1^2 + q3^2)
    R[2, 3] = 2 * (q2 * q3 + q1 * q0)
    R[3, 1] = 2 * (q1 * q3 + q2 * q0)
    R[3, 2] = 2 * (q2 * q3 - q1 * q0)
    R[3, 3] = 1 - 2 * (q1^2 + q2^2)
    return R
end

# Quaternion product
function quaternion_product(q::AbstractVector{T}, p::AbstractVector{T}) :: AbstractVector{T} where T <: Real
    q0, q1, q2, q3 = q
    p0, p1, p2, p3 = p
    s = q0 * p0 - q1 * p1 - q2 * p2 - q3 * p3
    x = q0 * p1 + q1 * p0 + q2 * p3 - q3 * p2
    y = q0 * p2 - q1 * p3 + q2 * p0 + q3 * p1
    z = q0 * p3 + q1 * p2 - q2 * p1 + q3 * p0
    return [s, x, y, z]
end


function state_jacobian6dof(x::AbstractVector{T}, u::AbstractVector{T}, params::Dict) :: AbstractMatrix{T} where T <: Real
    A = ForwardDiff.jacobian(x_state -> dynamics6dof(x_state, u, params), x)
    return A
end

function control_jacobian6dof(x::AbstractVector{T}, u::AbstractVector{T}, params::Dict) :: AbstractMatrix{T} where T <: Real
    B = ForwardDiff.jacobian(u_state -> dynamics6dof(x, u_state, params), u)
    return B
end

function initialize_trajectory6dof(params::Dict) :: Tuple{AbstractArray{Real,2}, AbstractArray{Real,2}}
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
        dx = dynamics6dof(xk, uk, params)
        x_init[k+1, :] = xk + params["dt"] * dx
        u_init[k, :] = uk
    end

    return x_init, u_init
end

# Helper function to compute skew-symmetric matrix of a vector TODO: move to utils
function skew_symmetric3d(v::AbstractVector{T}) :: AbstractMatrix{T} where T <: Real
    x, y, z = v
    S = zeros(T, 3, 3)
    S[1, 2] = -z
    S[1, 3] = y
    S[2, 1] = z
    S[2, 3] = -x
    S[3, 1] = -y
    S[3, 2] = x
    return S
end

function skew_symmetric4d(v::AbstractVector{T}) :: AbstractMatrix{T} where T <: Real
    x, y, z = v
    S = zeros(T, 4, 4)
    S[1, 2] = -x
    S[1, 3] = -y
    S[1, 4] = -z
    S[2, 1] = x
    S[2, 3] = z
    S[2, 4] = -y
    S[3, 1] = y
    S[3, 2] = -z
    S[3, 4] = x
    S[4, 1] = z
    S[4, 2] = y
    S[4, 3] = -x
    return S
end


end # module RocketDynamics6dof