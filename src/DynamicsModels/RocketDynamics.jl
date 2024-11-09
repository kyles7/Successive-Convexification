module RocketDynamics

using ..AbstractDynamicsModel

export RocketDynamics, dynamics, state_jacobian, control_jacobian, initialize_trajectory

struct RocketDynamics <: AbstractDynamicsModel.DynamicsModel
    params::Dict
end

"""
    dynamics(model::RocketDynamics, x::Vector, u::Vector, params::Dict) -> Vector

Compute the dynamics of the rocket.

# Returns
- `dx::Vector`: Derivative of the state vector.
"""
function dynamics(model::RocketDynamics, x::Vector, u::Vector, params::Dict)
    # Extract state variables
    px, py, vx, vy, m = x
    # Extract control inputs
    Tx, Ty = u
    # Parameters
    g = params["gravity"]
    
    # Compute dynamics
    dp = [vx, vy]                      # Position derivatives
    dv = [Tx / m, Ty / m - g]          # Velocity derivatives
    dm = -sqrt(Tx^2 + Ty^2) / params["I_sp"]  # Mass derivative
    
    dx = [dp; dv; dm]
    return dx
end

"""
    state_jacobian(model::RocketDynamics, x::Vector, u::Vector, params::Dict) -> Matrix

Compute the state Jacobian ∂f/∂x for the rocket dynamics.
"""
function state_jacobian(model::RocketDynamics, x::Vector, u::Vector, params::Dict)
    # Initialize Jacobian matrix
    n_states = length(x)
    A = zeros(n_states, n_states)
    px, py, vx, vy, m = x
    Tx, Ty = u
    g = params["gravity"]
    
    # ∂(dp_x)/∂(v_x)
    A[1, 3] = 1.0
    # ∂(dp_y)/∂(v_y)
    A[2, 4] = 1.0
    # ∂(dv_x)/∂(m)
    A[3, 5] = -Tx / m^2
    # ∂(dv_y)/∂(m)
    A[4, 5] = -Ty / m^2
    # ∂(dm)/∂(Tx) and ∂(dm)/∂(Ty) are zero with respect to state variables
    
    return A
end

"""
    control_jacobian(model::RocketDynamics, x::Vector, u::Vector, params::Dict) -> Matrix

Compute the control Jacobian ∂f/∂u for the rocket dynamics.
"""
function control_jacobian(model::RocketDynamics, x::Vector, u::Vector, params::Dict)
    n_states = length(x)
    n_controls = length(u)
    B = zeros(n_states, n_controls)
    px, py, vx, vy, m = x
    Tx, Ty = u
    g = params["gravity"]
    I_sp = params["I_sp"]
    
    # ∂(dv_x)/∂(T_x)
    B[3, 1] = 1 / m
    # ∂(dv_y)/∂(T_y)
    B[4, 2] = 1 / m
    
    # ∂(dm)/∂(T_x)
    T_norm = sqrt(Tx^2 + Ty^2)
    if T_norm != 0
        B[5, 1] = -Tx / (I_sp * T_norm)
        B[5, 2] = -Ty / (I_sp * T_norm)
    else
        B[5, 1] = 0.0
        B[5, 2] = 0.0
    end
    
    return B
end

"""
    initialize_trajectory(model::RocketDynamics, params::Dict) -> (Array, Array)

Initialize the state and control trajectories for the rocket dynamics model.
"""
function initialize_trajectory(model::RocketDynamics, params::Dict)
    N = params["N"]
    n_states = params["n_states"]
    n_controls = params["n_controls"]
    
    # Initialize state and control trajectories
    x_init = zeros(N, n_states)
    u_init = zeros(N-1, n_controls)
    
    # Set initial state
    x_init[1, :] = params["x0"]
    
    # Simple initial guess: zero control inputs and propagate dynamics
    for k in 1:N-1
        xk = x_init[k, :]
        uk = params["u_guess"]
        x_init[k+1, :] = xk + params["dt"] * dynamics(model, xk, uk, params)
        u_init[k, :] = uk
    end
    
    return x_init, u_init
end

end # module RocketDynamics
