# test/test_discretization.jl
module TestDynamics

export runTestDynamics

using Test
using YAML
using LinearAlgebra
using ..MainModule
using DifferentialEquations
using Plots

function runTestDynamics()
    config_path = joinpath(@__DIR__, "..", "..", "configs", "config6dof.yaml")
    params = Parameters.load_parameters(config_path)

    # normalize the parameters
    nondimensionalize!(params)

    # initialize the trajectory with the python library trajectory for K = 5
    x_ref, u_ref = initialize_trajectory6dof(params)
    x_ref = Array{Float64,2}(x_ref)
    u_ref = Array{Float64,2}(u_ref)
    # x_ref = [1.0 0.9466666666666668 0.8933333333333333 0.84 0.7866666666666666; 0.34815531191139565 0.27852424952911653 0.20889318714683738 0.13926212476455826 0.06963106238227913; 0.34815531191139565 0.27852424952911653 0.20889318714683738 0.13926212476455826 0.06963106238227913; 0.8703882797784891 0.6963106238227913 0.5222329678670934 0.34815531191139565 0.17407765595569782; -0.08703882797784891 -0.06963106238227913 -0.052223296786709346 -0.034815531191139566 -0.017407765595569783; -0.08703882797784891 -0.06963106238227913 -0.052223296786709346 -0.034815531191139566 -0.017407765595569783; -0.17407765595569782 -0.14100290132411525 -0.10792814669253264 -0.07485339206095007 -0.04177863742936748; 1.0 1.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
    # u_ref = [0 0 0 0 0; 0 0 0 0 0; 0.01360284592524594 0.01360284592524594 0.01360284592524594 0.01360284592524594 0.01360284592524594]
    
    println("x_1 = ", x_ref[:,1])
    println("u_1 = ", u_ref[:,1])
    f = dynamics6dof(x_ref[:,1], u_ref[:,1], params)

end # function runTestDiscretization

end # module TestDiscretization