# test/test_nomdim.jl

module TestNondimensionalization

export runTestNondimensionalization

using Test
using YAML
using LinearAlgebra
using ..MainModule
using DifferentialEquations
using Plots

function runTestNondimensionalization()
    config_path = joinpath(@__DIR__, "..", "..", "configs", "config6dof.yaml")
    params = Parameters.load_parameters(config_path)
    x_0, u_0 = initialize_trajectory6dof(params)
    x_0 = Array{Float64,2}(x_0)
    u_0 = Array{Float64,2}(u_0)
    println("Results of Dimensional initial trajectories: ")
    println("---------------------------------")
    println(x_0)
    println("---------------------------------")
    println(u_0)
    println("---------------------------------")
    
    nondimensionalize!(params)
    x_1, u_1 = initialize_trajectory6dof(params)
    x_1 = Array{Float64,2}(x_1)
    u_1 = Array{Float64,2}(u_1)
    # print the results
    # println("Results before redimensionalization: ")
    # println("---------------------------------")
    # println(x_1)
    # println("---------------------------------")
    # println(u_1)
    # println("---------------------------------")

    redimensionalize!(params)
    r_scale = norm(params["x0"][1:3])
    m_scale = params["m_wet"]

    for k in 1:size(x_1, 2)
        state_vec = @view x_1[:, k]
        x_redim!(state_vec, m_scale, r_scale)
        control_vec = @view u_1[:, k]
        u_redim!(control_vec, m_scale, r_scale)
    end

    # print the results
    # println("Results after redimensionalization: ")
    # println("---------------------------------")
    # println(x_1)
    # println("---------------------------------")
    # println(u_1)
    # println("---------------------------------")

    # Add a test compating 
    @test isapprox(norm(x_0), norm(x_1); rtol=1e-6)
    @test isapprox(norm(u_0), norm(u_1); rtol=1e-6)
    @test all(isapprox.(x_0, x_1; rtol=1e-6))
    @test all(isapprox.(u_0, u_1; rtol=1e-6))
end
end # module TestNondimensionalization
