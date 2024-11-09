# src/MainModule.jl
module MainModule
    include("Utilities/Parameters.jl")
    # include("DynamicsModels/AbstractDynamicsModel.jl")
    include("DynamicsModels/RocketDynamics.jl")

    using .Parameters, .RocketDynamics
    # using .AbstractDynamicsModel
    
    # Export everything from submodules
    export Parameters, RocketDynamics3dof, dynamics, state_jacobian, control_jacobian, initialize_trajectory
end
