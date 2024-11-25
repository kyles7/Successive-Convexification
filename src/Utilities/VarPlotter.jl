module VarPlotter

export plot_thrust_magnitude

using Plots
using LinearAlgebra
using CairoMakie

# function plot_thrust_magnitude(U::Array{T, 2}, params::Dict) where T <: Real
#     """
#     Plots the thrust magnitude for a given 3xN matrix of thrust vectors.

#     Parameters:
#     - thrust_vectors: 3xN matrix where each column is a thrust vector.
#     - T_max: Maximum thrust limit (horizontal line).
#     - T_min: Minimum thrust limit (horizontal line).
#     """
#     # Plot the thrust magnitude over time
#     thrust_magnitude = [norm(U[:, i]) for i in 1:size(U, 2)]
#     # Define the time steps
#     N = size(U, 2)
#     time_steps = 1:N
#     # Create the plot
#     plot(time_steps, thrust_magnitude, label="Thrust Magnitude", xlabel="Time Step", ylabel="Thrust Magnitude", linewidth=2)
#     println("T_min: ", params["T_min"])
#     println("T_max: ", params["T_max"])
#     plot!(time_steps, fill(params["T_max"], N), label="Max Thrust", linestyle=:dash, color=:red, linewidth=2)
#     plot!(time_steps, fill(params["T_min"], N), label="Min Thrust", linestyle=:dash, color=:green, linewidth=2)
#     display(plot)
#     # gui()
#     # show()
#     return nothing
# end



function plot_thrust_magnitude(U::Matrix{Float64}, params::Dict)
    # Compute thrust magnitude
    thrust_magnitude = [norm(U[:, i]) for i in 1:size(U, 2)]
    N = size(U, 2)
    time_steps = 1:N

    # Create a figure and axis
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Time Step", ylabel = "Thrust Magnitude")

    # Plot thrust magnitude
    lines!(ax, time_steps, thrust_magnitude, label = "Thrust Magnitude", color = :blue, linewidth = 2)

    # Plot max and min thrust as horizontal lines
    hlines!(ax, [params["T_max"]], label = "Max Thrust", color = :red, linestyle = :dash)
    hlines!(ax, [params["T_min"]], label = "Min Thrust", color = :green, linestyle = :dash)

    # Add legend and display the plot
    axislegend(ax)
    display(fig)
end



end  # module VarPlotter