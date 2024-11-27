module TrajectoryVisualization

export plot_trajectory

using LinearAlgebra
using CairoMakie
using GeometryBasics
using GLMakie
using ..MainModule


function plot_trajectory(X::Array{Float64,2}, U::Array{Float64,2}, thrust_scale::Float64, attitude_scale)
    # X ∈ ℝ^(14 x N) is the state trajectory
    # U ∈ ℝ^(3 x N) is the control trajectory

    fig = Figure(size = (800, 600))
    ax = Axis3(fig[1, 1], xlabel = "X, east", ylabel = "Y, north", zlabel = "Z, up", aspect = :data)

    # Extract positions over time
    rx = X[2, :]
    ry = X[3, :]
    rz = X[4, :]

    # Plot the trajectory path
    lines!(ax, rx, ry, rz, color = :black, linewidth = 1.5, linestyle = :dash, label = "Trajectory")

    # Loop over time steps to plot attitude and thrust vectors
    for k in 1:size(X, 2)  # size(X', 1) is N, the number of time steps
        # Extract quaternion components
        qw, qx, qy, qz = X[8, k], X[9, k], X[10, k], X[11, k]

        # Compute rotation matrix
        CBI = quaternion_to_rotation_matrix(qw, qx, qy, qz)

        # Compute attitude and thrust vectors
        attitude_vector = CBI' * [0.0, 0.0, 1.0]
        thrust_vector = CBI' * U[:, k]

        # Scale the vectors
        attitude_vector_scaled = attitude_vector .* attitude_scale
        thrust_vector_scaled = -thrust_vector .* thrust_scale

        # Convert to Vec3f0
        attitude_vec3 = Vec3f0(attitude_vector_scaled...)
        thrust_vec3 = Vec3f0(thrust_vector_scaled...)

        # Plot vectors at every Nth time step to reduce clutter
        if k % 2 == 0  # Adjust the step as needed
            position = Point3f0(rx[k], ry[k], rz[k])
            scatter!(ax, position, color = :blue)
            # Plot attitude vector
            arrows!(ax, [position], [attitude_vec3], color = :blue, linewidth = 1.3)
            # Plot thrust vector
            arrows!(ax, [position], [thrust_vec3], color = :red, linewidth = 1.3)
        end
    end
    
    # Optional: Plot the initial and final positions
    scatter!(ax, [rx[1]], [ry[1]], [rz[1]], color = :green, markersize = 10, label = "Start")
    scatter!(ax, [rx[end]], [ry[end]], [rz[end]], color = :red, markersize = 10, label = "End")

    # Create a disk at z = 0
    disk = HyperSphere(Point3f0(0, 0, 0), 20.0f0)

    # Plot the disk
    mesh!(ax, disk, color = (:lightgray, 0.5), transparency = true)

    # Set the plot title
    ax.title = "Trajectory Visualization"

    # Optionally adjust axis limits
    ax.limits = ((minimum(rx) - 50, maximum(rx) + 50),
                 (minimum(ry) - 50, maximum(ry) + 50),
                 (minimum(rz) - 50, maximum(rz) + 50))

    axislegend(ax)
    display(fig)
end # function plot_trajectory

end  # module TrajectoryVisualization
