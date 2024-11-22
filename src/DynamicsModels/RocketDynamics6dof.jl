module RocketDynamics6dof

using LinearAlgebra
using ForwardDiff
using DifferentialEquations

export RocketDynamics_6dof, dynamics6dof, state_jacobian6dof, control_jacobian6dof, initialize_trajectory6dof, quaternion_to_rotation_matrix, calculate_discretization, x_nondim!, u_nondim!, x_redim!, u_redim!, nondimensionalize!, redimensionalize!

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
    rTB = params["r_T_B"]
    # Extract control inputs
    TBx, TBy, TBz = u[1:3]
    # Compute thrust and moment vectors
    FIx, FIy, FIz = quaternion_to_rotation_matrix(q0, q1, q2, q3) * [TBx; TBy; TBz]
    MBx, MBy, MBz = skew_symmetric3d(rTB) * [TBx; TBy; TBz]

    # # Thrust force in body frame
    # Fbx, Fby, Fbz = u[1:3]
    # # Thrust moment in body frame
    # Mbx, Mby, Mbz = u[4:6]

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
    return R'
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
    """
    Initialize the trajectory.

    :param params: Dictionary containing parameters like T_max and T_min
    :return: The initialized X (n_states x K) and U (n_controls x K)
    """
    N = params["N"]
    n_states = params["n_states"]
    n_controls = params["n_controls"]
    T_max = params["T_max"]
    T_min = params["T_min"]

    # Initialize state and control matrices
    X = zeros(n_states, N) # 14 x 10
    U = zeros(n_controls, N) # 3 x 10 

    m0 = params["m_wet"]
    r0 = params["x0"][1:3]
    v0 = params["x0"][4:6]
    q0 = params["x0"][7:10]
    omega0 = params["x0"][11:13]
    mf = params["m_dry"]
    rf = params["xf"][1:3]
    vf = params["xf"][4:6]
    qf = params["xf"][7:10]
    omegaf = params["xf"][11:13]

    for k in 1:N
        alpha1 = (N-k) / N
        alpha2 = k / N

        m_k = alpha1 * m0 + alpha2 * mf
        r_I_k = alpha1 * r0 + alpha2 * rf
        v_I_k = alpha1 * v0 + alpha2 * vf
        q_B_I_k = alpha1 * q0 + alpha2 * qf
        omega_B_k = alpha1 * omega0 + alpha2 * omegaf

        X[:, k] = vcat(r_I_k, v_I_k, q_B_I_k, omega_B_k, m_k)
        U[:, k] = ((T_max - T_min) /2) * [0.0, 0.0, 1.0]
    end

    return X, U
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

function calculate_discretization(X::AbstractMatrix{T}, U::AbstractMatrix{T}, sigma, params::Dict) :: Tuple{AbstractMatrix{T}, AbstractMatrix{T}, AbstractMatrix{T}, AbstractMatrix{T}, AbstractMatrix{T}} where T <: Real
   """
    Calculate discretization for given states, inputs, and total time.

    :param X: Matrix of states for all time points (n_x x K)
    :param U: Matrix of inputs for all time points (n_u x K)
    :param sigma: Total time
    :param params: Parameters dictionary
    :return: Discretization matrices A_bar, B_bar, C_bar, S_bar, z_bar
    """
    # Extract dimensions
    
    K = size(X, 2)
    n_x = size(X, 1)
    n_u = size(U, 1)
    dt = sigma / (K - 1)

    #preallocate matrices
    A_bar = zeros(n_x * n_x, K-1)
    B_bar = zeros(n_x * n_u, K-1)
    C_bar = zeros(n_x * n_u, K-1)
    S_bar = zeros(n_x, K-1)
    z_bar = zeros(n_x, K-1)

    #indices for augmented state vec V
    x_index = 1:n_x
    A_bar_indices = maximum(x_index) + 1:maximum(x_index) + n_x * n_x #??
    B_bar_indices = maximum(A_bar_indices) + 1:maximum(A_bar_indices) + n_x * n_u
    C_bar_indices = maximum(B_bar_indices) + 1:maximum(B_bar_indices) + n_x * n_u
    S_bar_indices = maximum(C_bar_indices) + 1:maximum(C_bar_indices) + n_x
    z_bar_indices = maximum(S_bar_indices) + 1:maximum(S_bar_indices) + n_x

    V0_length = maximum(z_bar_indices)

    #preallocate augmented state vec V
    V0 = zeros(V0_length)
    V0[A_bar_indices] = vec(Matrix{Float64}(I, n_x, n_x))

    for k in 1:K-1
        #set initial augmented state
        V0[x_index] = X[:, k]

        # define ODE function
        # inputs vector to store derivatives, current state vec, params, current time
        # function f!(dVdt, V, p, t)
        #     # extract variables from V
        #     xk = V[x_index]
        #     Phi = reshape(V[A_bar_indices], n_x, n_x)
        #     B_bar_V = reshape(V[B_bar_indices], n_x, n_u)
        #     C_bar_V = reshape(V[C_bar_indices], n_x, n_u)
        #     S_bar_V = V[S_bar_indices]
        #     z_bar_V = V[z_bar_indices]

        #     # interpolate control input
        #     alpha = t / dt
        #     u = (1 - alpha) * U[:, k] + alpha * U[:, k+1]

        #     # compute dynamics and jacobian
        #     f = dynamics6dof(xk, u, params)
        #     A = state_jacobian6dof(xk, u, params)
        #     B = control_jacobian6dof(xk, u, params)

        #     # compute derivatives
        #     dxdt = f  
        #     dPhidt = A * Phi
        #     dB_bar_dt = A * B_bar_V + B
        #     dC_bar_dt = A * C_bar_V
        #     dS_bar_dt = A * S_bar_V
        #     dz_bar_dt = A * z_bar_V + f

        #     # pack derivatives into dVdt
        #     dVdt = zeros(V0_length)
        #     dVdt[x_index] = dxdt
        #     dVdt[A_bar_indices] = vec(dPhidt)
        #     dVdt[B_bar_indices] = vec(dB_bar_dt)
        #     dVdt[C_bar_indices] = vec(dC_bar_dt)
        #     dVdt[S_bar_indices] = dS_bar_dt
        #     dVdt[z_bar_indices] = dz_bar_dt

        #     return nothing
        # end
        function f!(dVdt, V, sigma, t)
            alpha = (dt - t) / dt
            beta = t / dt
            x = V[x_index]
            u = U[:, k] + (t / dt) * (U[:, k+1] - U[:, k])
            Phi_A_xi = inv(reshape(V[A_bar_indices], n_x, n_x))
            A_subs = sigma * state_jacobian6dof(x, u, params)
            B_subs = sigma * control_jacobian6dof(x, u, params)
            f_subs = dynamics6dof(x, u, params)
            # initialize dVdt
            dVdt = zeros(V0_length)
            dVdt[x_index] = sigma * f_subs'
            dVdt[A_bar_indices] = vec(A_subs * reshape(V[A_bar_indices], n_x, n_x))
            dVdt[B_bar_indices] = vec(Phi_A_xi * B_subs) * alpha
            dVdt[C_bar_indices] = vec(Phi_A_xi * B_subs) * beta
            dVdt[S_bar_indices] = (Phi_A_xi * f_subs)'
            z_t = -A_subs * x - B_subs * u
            dVdt[z_bar_indices] = (Phi_A_xi * z_t)
            return dVdt
        end

        #test = f!(zeros(V0_length), V0, sigma, 0.0) # using this allows you to step thru f!
        # Integrate ODE from t=0 to t=dt
        tspan = (0.0, dt)
        prob = ODEProblem(f!, V0, tspan, sigma)
        sol = solve(prob, AutoVern7(Rodas5()), reltol = 1e-8, abstol = 1e-8; saveat=dt)
        #get V at t = dt
        V_end = sol.u[end]

        # Extract Phi and other matrices
        Phi = reshape(V_end[A_bar_indices], n_x, n_x)
        A_bar[:, k] = vec(Phi)

        B_bar_V = reshape(V_end[B_bar_indices], n_x, n_u)
        B_bar[:, k] = vec(Phi * B_bar_V)

        C_bar_V = reshape(V_end[C_bar_indices], n_x, n_u)
        C_bar[:, k] = vec(Phi * C_bar_V)

        S_bar_V = V_end[S_bar_indices]
        S_bar[:, k] = Phi * S_bar_V

        z_bar_V = V_end[z_bar_indices]
        z_bar[:, k] = Phi * z_bar_V

    end
    return A_bar, B_bar, C_bar, S_bar, z_bar
end

function x_nondim!(x::AbstractVector{T}, m_scale::T, r_scale::T) :: Nothing where T <: Real
    """
        Nondimensionalize the state vector x.
    """
    # Nondimensionalize position
    x[1:3] /= r_scale
    # Nondimensionalize velocity
    x[4:6] /= r_scale
    # Nondimensionalize mass
    x[14] /= m_scale
    return nothing
end

function u_nondim!(u::T, m_scale::T, r_scale::T) :: T where T <: Real
    """
    Nondimensionalize the control vector u.
    """
    u /= (m_scale * r_scale)
    return u
end

function x_redim!(x::AbstractVector{T}, m_scale::T, r_scale::T) :: Nothing where T <: Real
    """
    Redimensionalize the state vector x.
    """
    # Redimensionalize position
    x[1:3] *= r_scale
    # Redimensionalize velocity
    x[4:6] *= r_scale
    return nothing
end

function u_redim!(u::T, m_scale::T, r_scale::T) :: T where T <: Real
    """
    Redimensionalize the control vector u.
    """
    u *= (m_scale * r_scale)
    return u
end

# Nondimensionalization and redimensionalization functions
function nondimensionalize!(params::Dict)
    """
    Nondimensionalize the parameters in the params dictionary.
    """
    # Nondimensionalize parameters
    r_scale = norm(params["x0"][1:3])
    m_scale = params["m_wet"]

    params["r_T_B"] /= r_scale
    params["gravity_vector"] /= r_scale
    params["inertia_matrix"] /= m_scale * r_scale^2

    # Nondimensionalize initial and final states
    x_nondim!(params["x0"], m_scale, r_scale)
    x_nondim!(params["xf"], m_scale, r_scale)

    params["T_max"] = u_nondim!(params["T_max"], m_scale, r_scale)
    params["T_min"] = u_nondim!(params["T_min"], m_scale, r_scale)

    # Nondimensionalize masses
    params["m_wet"] /= m_scale
    params["m_dry"] /= m_scale
    return nothing
end

function redimensionalize!(params::Dict)
    """
    Redimensionalize the parameters in the params dictionary.
    """
    # Redimensionalize parameters
    r_scale = norm(params["x0"][1:3])
    m_scale = params["m_wet"]

    params["r_T_B"] *= r_scale
    params["gravity_vector"] *= r_scale
    params["inertia_matrix"] *= m_scale * r_scale^2

    # Redimensionalize initial and final states
    x_redim(params["x0"])
    x_redim(params["xf"])

    # Redimensionalize control limits
    u_redim(params["T_max"])
    u_redim(params["T_min"])

    # Redimensionalize masses
    params["m_wet"] *= m_scale
    params["m_dry"] *= m_scale
    return nothing
end



end # module RocketDynamics6dof