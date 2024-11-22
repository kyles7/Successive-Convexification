# test/test_discretization.jl
module TestDiscretization

export runTestDiscretization

using Test
using YAML
using LinearAlgebra
using ..MainModule
using DifferentialEquations
using Plots

function runTestDiscretization()
    config_path = joinpath(@__DIR__, "..", "..", "configs", "config6dof.yaml")
    params = Parameters.load_parameters(config_path)

    # normalize the parameters
    nondimensionalize!(params)

    x_ref, u_ref = initialize_trajectory6dof(params)
    x_ref = Array{Float64,2}(x_ref)
    u_ref = Array{Float64,2}(u_ref)

    # println("x_ref: ", x_ref)
    # println("u_ref: ", u_ref)

    # make discretization matrices
    A_bar, B_bar, C_bar, S_bar, z_bar = calculate_discretization(x_ref, u_ref, params["sigma"], params)
    
    println("A_bar: ", A_bar)
    println("B_bar: ", B_bar)
    println("C_bar: ", C_bar)
    println("S_bar: ", S_bar)
    println("z_bar: ", z_bar)


end # function runTestDiscretization

end # module TestDiscretization