# src/DynamicsModels/AbstractDynamicsModel.jl

module AbstractDynamicsModel

export DynamicsModel, dynamics, state_jacobian, control_jacobian, initialize_trajectory

"""
    abstract type DynamicsModel

An abstract type representing a general dynamics model. All specific dynamics models
(e.g., RocketDynamics, QuadrotorDynamics) should subtype this abstract type.
"""
abstract type DynamicsModel end

"""
    dynamics(model::DynamicsModel, x::Vector, u::Vector, params::Dict) -> Vector

Compute the dynamics (state derivatives) given the current state `x`, control input `u`,
and system parameters `params`.

# Arguments
- `model::DynamicsModel`: An instance of a dynamics model.
- `x::Vector`: State vector at the current time step.
- `u::Vector`: Control input vector at the current time step.
- `params::Dict`: Dictionary of system parameters.

# Returns
- `dx::Vector`: Derivative of the state vector.
"""
function dynamics(model::DynamicsModel, x::Vector, u::Vector, params::Dict)
    throw(NotImplementedError("The `dynamics` function must be implemented for $(typeof(model))"))
end

"""
    state_jacobian(model::DynamicsModel, x::Vector, u::Vector, params::Dict) -> Matrix

Compute the Jacobian of the dynamics with respect to the state vector `x`.

# Arguments
- `model::DynamicsModel`: An instance of a dynamics model.
- `x::Vector`: State vector at the current time step.
- `u::Vector`: Control input vector at the current time step.
- `params::Dict`: Dictionary of system parameters.

# Returns
- `A::Matrix`: State Jacobian matrix ∂f/∂x.
"""
function state_jacobian(model::DynamicsModel, x::Vector, u::Vector, params::Dict)
    throw(NotImplementedError("The `state_jacobian` function must be implemented for $(typeof(model))"))
end

"""
    control_jacobian(model::DynamicsModel, x::Vector, u::Vector, params::Dict) -> Matrix

Compute the Jacobian of the dynamics with respect to the control vector `u`.

# Arguments
- `model::DynamicsModel`: An instance of a dynamics model.
- `x::Vector`: State vector at the current time step.
- `u::Vector`: Control input vector at the current time step.
- `params::Dict`: Dictionary of system parameters.

# Returns
- `B::Matrix`: Control Jacobian matrix ∂f/∂u.
"""
function control_jacobian(model::DynamicsModel, x::Vector, u::Vector, params::Dict)
    throw(NotImplementedError("The `control_jacobian` function must be implemented for $(typeof(model))"))
end

"""
    initialize_trajectory(model::DynamicsModel, params::Dict) -> (Array, Array)

Initialize the state and control trajectories for the optimization algorithm.

# Arguments
- `model::DynamicsModel`: An instance of a dynamics model.
- `params::Dict`: Dictionary of system parameters.

# Returns
- `(x_init::Array, u_init::Array)`: Tuple containing the initial state and control trajectories.
"""
function initialize_trajectory(model::DynamicsModel, params::Dict)
    throw(NotImplementedError("The `initialize_trajectory` function must be implemented for $(typeof(model))"))
end

end # module AbstractDynamicsModel
