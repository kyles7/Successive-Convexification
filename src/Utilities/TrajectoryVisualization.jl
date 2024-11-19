using LinearAlgebra
using CairoMakie
using GeometryBasics
using GLMakie

# TODO: Make this a module and export the function

# TODO: get this function from the RocketDynamics6dof.jl file
function quaternion_to_rotation_matrix(qw, qx, qy, qz)
    return [
        1 - 2 * (qy^2 + qz^2)     2 * (qx * qy + qw * qz)     2 * (qx * qz - qw * qy);
        2 * (qx * qy - qw * qz)   1 - 2 * (qx^2 + qz^2)       2 * (qy * qz + qw * qx);
        2 * (qx * qz + qw * qy)   2 * (qy * qz - qw * qx)     1 - 2 * (qx^2 + qy^2)
    ]
end

function plot_trajectory(X::Array{Float64,2}, U::Array{Float64,2}, thrust_scale::Float64, attitude_scale)
    fig = Figure(size = (800, 600))
    ax = Axis3(fig[1, 1], xlabel = "X, east", ylabel = "Y, north", zlabel = "Z, up")

    # Extract positions over time
    rx = X[:, 1]
    ry = X[:, 2]
    rz = X[:, 3]

    # Plot the trajectory path
    lines!(ax, rx, ry, rz, color = :black, linewidth = 1.5, linestyle = :dash, label = "Trajectory")


    # Loop over time steps to plot attitude and thrust vectors
    for i in 1:size(X, 1)  # size(X, 1) is N, the number of time steps
        # Extract quaternion components
        qw, qx, qy, qz = X[i, 7], X[i, 8], X[i, 9], X[i, 10]

        # Compute rotation matrix
        CBI = quaternion_to_rotation_matrix(qw, qx, qy, qz)

        # Compute attitude and thrust vectors
        attitude_vector = CBI' * [0.0, 0.0, 1.0]
        thrust_vector = CBI' * U[i, :]

        # Scale the vectors
        attitude_vector_scaled = attitude_vector .* attitude_scale
        thrust_vector_scaled = -thrust_vector .* thrust_scale

        # Convert to Vec3f0
        attitude_vec3 = Vec3f0(attitude_vector_scaled...)
        thrust_vec3 = Vec3f0(thrust_vector_scaled...)

        # Plot vectors at every Nth time step to reduce clutter
        if i % 1 == 0  # Adjust the step as needed
            position = Point3f0(rx[i], ry[i], rz[i])
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

    # Add legend
    axislegend(ax)

    # Display the figure
    display(fig)
end


# Main
# Example data with N time steps
N = 10

# Constants for scaling
thrust_scale = 0.0002
attitude_scale = 10

# Generate example trajectory data (replace with your actual data)
X = zeros(N, 14)
U = zeros(N, 3)
X = [0.0 200.0 500.0 -50.0 -100.0 -50.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 30000.0;
     -5.0 190.0 495.0 -50.0 -100.0 -50.817165 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29998.33276399178;
     -10.0 180.0 489.9182835 -50.0 -100.0 -51.63432091305876 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29996.66552798356;
    -15.0 170.0 484.7548514086941 -50.0 -100.0 -52.451467738166144 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29994.998291975342;
    -20.0 160.0 479.5097046348775 -50.0 -100.0 -53.26860547431188 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29993.331055967123;
    -25.0 150.0 474.1828440874463 -50.0 -100.0 -54.08573412048552 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29991.663819958903;
    -30.0 140.0 468.77427067539776 -50.0 -100.0 -54.902853675676425 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29989.996583950684;
    -35.0 130.0 463.2839853078301 -50.0 -100.0 -55.71996413887381 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29988.329347942465;
    -40.0 120.0 457.7119888939427 -50.0 -100.0 -56.53706550906672 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29986.662111934245;
    -45.0 110.0 452.058282343036 -50.0 -100.0 -57.35415778524402 1.0 0.0 0.0 0.0 0.0 0.0 0.0 29984.994875926026]

U =  [0.0 0.0 49050.0;
        0.0 0.0 49050.0;
        0.0 0.0 49050.0;
        0.0 0.0 49050.0;
        0.0 0.0 49050.0;
        0.0 0.0 49050.0;
        0.0 0.0 49050.0;
        0.0 0.0 49050.0;
        0.0 0.0 49050.0;
        0.0 0.0 49050.0]


# Plot the trajectory
plot_trajectory(X, U, thrust_scale, attitude_scale)
