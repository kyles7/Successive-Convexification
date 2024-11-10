# src/MainModule.jl
module MainModule
    ######## use this block for test_dynamics ########
    # include("Utilities/Parameters.jl")  
    # include("DynamicsModels/RocketDynamics.jl")  # use this for test_dynamics
    # using .Parameters, .RocketDynamics
    # export Parameters, RocketDynamics3dof, dynamics, state_jacobian, control_jacobian, initialize_trajectory

    ######## use this block for test_6dofdynamics #######
    include("Utilities/Parameters.jl")
    include("DynamicsModels/RocketDynamics6dof.jl")  # use this for test_6dofdynamics
    using .Parameters, .RocketDynamics
    export Parameters, RocketDynamics_6dof, dynamics, state_jacobian, control_jacobian, initialize_trajectory

#   include("Utilities/Parameters.jl")  
#     # include("DynamicsModels/AbstractDynamicsModel.jl")
#     #include("DynamicsModels/RocketDynamics.jl")  # use this for test_dynamics
#     include("DynamicsModels/RocketDynamics6dof.jl")  # use this for test_6dofdynamics

#     using .Parameters, .RocketDynamics
#     # using .AbstractDynamicsModel
    
#     # Export everything from submodules
#     export Parameters, RocketDynamics3dof, dynamics, state_jacobian, control_jacobian, initialize_trajectory
end