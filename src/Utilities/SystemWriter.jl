module SystemWriter

export write_trajectory_data

using HDF5

function write_trajectory_data(file_path::String, X_all::Vector{Matrix{Float64}}, U_all::Vector{Matrix{Float64}}, sigma_all::Vector{Float64}, params::Dict)
    h5open(file_path, "w") do file
        # Create groups for state and control trajectories
        group_X = create_group(file, "X_all")
        group_U = create_group(file, "U_all")
        group_sigma = create_group(file, "sigma_all")
        group_params = create_group(file, "params")

        # Write each trajectory to the respective group
        for (i, X) in enumerate(X_all)
            group_X["iteration_$i"] = X  # Store 14xK matrix for states
        end
        for (i, U) in enumerate(U_all)
            group_U["iteration_$i"] = U  # Store 3xK matrix for controls
        end

        # Write sigma values
        for (i, sigma) in enumerate(sigma_all)
            group_sigma["iteration_$i"] = sigma
        end

    end
end

end 


