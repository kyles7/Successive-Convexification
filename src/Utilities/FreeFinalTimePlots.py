import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py
import matplotlib.lines as mlines

def quaternion_to_rotation_matrix(q0, q1, q2, q3):
    R = np.zeros((3, 3))
    R[0, 0] = 1 - 2 * (q2**2 + q3**2)
    R[0, 1] = 2 * (q1 * q2 + q0 * q3)
    R[0, 2] = 2 * (q1 * q3 - q0 * q2)
    R[1, 0] = 2 * (q1 * q2 - q0 * q3)
    R[1, 1] = 1 - 2 * (q1**2 + q3**2)
    R[1, 2] = 2 * (q2 * q3 + q0 * q1)
    R[2, 0] = 2 * (q1 * q3 + q0 * q2)
    R[2, 1] = 2 * (q2 * q3 - q0 * q1)
    R[2, 2] = 1 - 2 * (q1**2 + q2**2)
    return R

def compute_tilt_angles(quaternions):
    """
    Computes tilt angles for a batch of quaternions in 4xN format.
    
    Parameters:
        quaternions: A 2D array of shape (4, N), where each column is [q0, q1, q2, q3].
        
    Returns:
        tilt_angles: A 1D array of tilt angles (in radians) for each quaternion.
    """
    # # Normalize quaternions
    # norms = np.linalg.norm(quaternions, axis=0, keepdims=True)  # Norm along columns
    # quaternions = quaternions / norms  # Broadcasting normalization

    # Extract quaternion components
    q0, q1, q2, q3 = quaternions  # Rows correspond to q0, q1, q2, q3

    # Compute the body-fixed Z-axis in the inertial frame for each quaternion
    z_body_in_inertial = np.array([
        2 * (q1 * q3 - q0 * q2),
        2 * (q2 * q3 + q0 * q1),
        1 - 2 * (q1**2 + q2**2)
    ])

    # Compute the dot product with the inertial Z-axis [0, 0, 1]
    cos_theta = z_body_in_inertial[2, :]  # Dot product is just the Z-component

    # Clamp the values to the range [-1, 1] for numerical stability
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    # Compute tilt angles
    tilt_angles = np.arccos(cos_theta)

    return tilt_angles


def read_trajectory_data(file_path):
    X_all = []
    U_all = []
    sigma_all = []
    with h5py.File(file_path, 'r') as file:
        # Read all state trajectories
        for key in file['X_all']:
            X_all.append(np.array(file['X_all'][key]).T)
        
        # Read all control trajectories
        for key in file['U_all']:
            U_all.append(np.array(file['U_all'][key]).T)
    
        # Read the sigma values
        for key in file['sigma_all']:
            sigma_all.append(np.array(file['sigma_all'][key]))
        # sigma = np.array(file['sigma_all'])

    return X_all, U_all, sigma_all

def plot_trajectory(X, U, iter_label, thrust_scale, attitude_scale, ax=None):
    """
    Plots a 3D trajectory with attitude and thrust vectors.

    Args:
        X (np.ndarray): State trajectory (14 x N).
        U (np.ndarray): Control trajectory (3 x N).
        iter_label (int): Iteration number for the plot title.
        thrust_scale (float): Scaling factor for thrust vectors.
        attitude_scale (float): Scaling factor for attitude vectors.
    """
    if ax is None:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("X, east")
    ax.set_ylabel("Y, north")
    ax.set_zlabel("Z, up")
    ax.set_title("Optimized Trajectory for Iteration ${}$".format(iter_label))

    # Extract positions over time
    rx = X[1, :]  # X[2, :] in Julia (1-based indexing)
    ry = X[2, :]  # X[3, :] in Julia
    rz = X[3, :]  # X[4, :] in Julia

    # Plot the trajectory path
    ax.plot(rx, ry, rz, color='black', linestyle='--', linewidth=1, label="Trajectory")

    # Loop over time steps to plot attitude and thrust vectors
    for k in range(X.shape[1]):
        # Extract quaternion components
        qw, qx, qy, qz = X[7, k], X[8, k], X[9, k], X[10, k]  # X[8:11, k] in Julia

        # Compute rotation matrix
        CBI = quaternion_to_rotation_matrix(qw, qx, qy, qz)

        # Compute attitude and thrust vectors
        attitude_vector = CBI.T @ np.array([0.0, 0.0, 1.0])
        thrust_vector = CBI.T @ U[:, k]

        # Scale the vectors
        attitude_vector_scaled = attitude_vector * attitude_scale
        thrust_vector_scaled = -thrust_vector * thrust_scale

        # Plot vectors at every nth time step to reduce clutter
        if k % 1 == 0:
            position = np.array([rx[k], ry[k], rz[k]])
            ax.scatter(*position, color='blue', s=5)

            # Plot attitude vector
            ax.plot([position[0], position[0] + attitude_vector_scaled[0]],
                [position[1], position[1] + attitude_vector_scaled[1]],
                [position[2], position[2] + attitude_vector_scaled[2]], color='blue', linewidth=1.3)

            # Plot thrust vector
            ax.plot([position[0], position[0] + thrust_vector_scaled[0]],
                [position[1], position[1] + thrust_vector_scaled[1]],
                [position[2], position[2] + thrust_vector_scaled[2]], color='red', linewidth=1.3)

    # Optional: Plot the initial and final positions
    ax.scatter(rx[0], ry[0], rz[0], color='green', s=10, label="Start")
    ax.scatter(rx[-1], ry[-1], rz[-1], color='red', s=10, label="End")

    # Adjust axis limits
    # Equalize the axes
    max_range = max(max(rx) - min(rx), max(ry) - min(ry), max(rz) - min(rz)) / 2.0
    mid_x = (max(rx) + min(rx)) / 2.0
    mid_y = (max(ry) + min(ry)) / 2.0
    mid_z = (max(rz) + min(rz)) / 2.0

    # ax.set_xlim(mid_x - max_range, mid_x + max_range)
    # ax.set_ylim(mid_y - max_range, mid_y + max_range)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)
    # ax.set_xlim([min(rx) - 50, max(rx) + 50])
    # ax.set_ylim([min(ry) - 50, max(ry) + 50])
    # ax.set_zlim([min(rz) - 50, max(rz) + 50])

    # Show legend and plot
    ax.set_aspect('equal')
    ax.legend()
    # plt.show()

def plot_trajectory_2D(X_all, U_all, iter_index, thrust_scale, attitude_scale, pos_indices, vec_indices, ax=None):
    """
    Plots the 2D trajectory and vectors for any given viewpoint.
    
    Parameters:
    - X_all: List of state trajectories (each element is a matrix of shape 14xK).
    - U_all: List of control trajectories (each element is a matrix of shape 3xK).
    - iter_index: Index of the iteration to plot.
    - thrust_scale: Scaling factor for the thrust vector.
    - attitude_scale: Scaling factor for the attitude vector.
    - pos_indices: Tuple of two indices for the position (e.g., (1, 3) for East-Up).
    - vec_indices: Tuple of two indices for the vectors (e.g., (0, 2) for thrust/attitude in East-Up).
    """
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

    # Extract the indices for positions
    pos1, pos2 = pos_indices
    vec1, vec2 = vec_indices

    labels = ["East - Position (m)", "North - Position (m)", "Up - Position (m)"]
    # Plot the trajectory in the specified plane
    # ax.plot(X_all[iter_index][pos1, :], X_all[iter_index][pos2, :], label="Trajectory", color='black', linewidth=0.9, linestyle='--')
    ax.plot(X_all[iter_index][pos1, :], X_all[iter_index][pos2, :], label="Trajectory", color='black', marker='o', markersize=3, linestyle='--', linewidth=0.85)

    for k in range(X_all[iter_index].shape[1]):
        # Extract quaternion components
        qw, qx, qy, qz = X_all[iter_index][7, k], X_all[iter_index][8, k], X_all[iter_index][9, k], X_all[iter_index][10, k]

        # Compute rotation matrix
        CBI = quaternion_to_rotation_matrix(qw, qx, qy, qz)

        # Compute attitude and thrust vectors
        attitude_vector = CBI.T @ np.array([0.0, 0.0, 1.0])
        thrust_vector = CBI.T @ U_all[iter_index][:, k]

        # Scale the vectors
        attitude_vector_scaled = attitude_vector * attitude_scale
        thrust_vector_scaled = -thrust_vector * thrust_scale

        # Plot vectors at every nth time step to reduce clutter
        if k % 1 == 0:  # Adjust the step to reduce clutter
            position = np.array([X_all[iter_index][pos1, k], X_all[iter_index][pos2, k]])
            
            # Plot attitude vector
            ax.plot(
                [position[0], position[0] + attitude_vector_scaled[vec1]],
                [position[1], position[1] + attitude_vector_scaled[vec2]],
                color='blue', linewidth=1.5
            )

            # Plot thrust vector
            ax.plot(
                [position[0], position[0] + thrust_vector_scaled[vec1]],
                [position[1], position[1] + thrust_vector_scaled[vec2]],
                color='red', linewidth=1.5
            )
    
    # ax.set_xlabel(f"Position {pos1}")
    ax.set_xlabel(labels[pos1 - 1])
    ax.set_ylabel(labels[pos2 - 1])
    ax.set_title(f"2D Trajectory View for Iteration {iter_index + 1}")
    # ax.set_ylabel(f"Position {pos2}")
    # ax.set_title(f"Position {pos1} vs {pos2}")
    ax.grid(True)
    # Add a blue line with the label "Thrust" to the legend
    thrust_line = mlines.Line2D([], [], color='red', label='Thrust')
    ax.add_line(thrust_line)
    # Add a red line with the label "Attitude" to the legend
    attitude_line = mlines.Line2D([], [], color='blue', label='Attitude')
    ax.add_line(attitude_line)
    # ax.set_aspect('equal')
    ax.legend()


    rx = X_all[iter_index][pos1, :]
    ry = X_all[iter_index][pos2, :]
    max_range = max(max(rx) - min(rx), max(ry) - min(ry)) / 2.0
    mid_x = (max(rx) + min(rx)) / 2.0
    mid_y = (max(ry) + min(ry)) / 2.0

    ax.set_xlim(mid_x - max_range-20, mid_x + max_range+20)
    ax.set_ylim(mid_y - max_range-20, mid_y + max_range+40)
    # plt.show()



def plot_state_var(vec, label, units=None, lower_constraint=None, upper_constraint=None, lower_constraint_label=None, upper_constraint_label=None, sigma=15, ax=None):
    # create new figure
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
    # convert time steps to true time with the final time sigma
    time = np.linspace(0, sigma, len(vec))
    ax.plot(time, vec, label=label, color='blue', linewidth=1.3, marker='o')  # plot the data
    # ax.plot(vec, label=label, color='blue', linewidth=1.3, marker='o')  # plot the data
    ax.set_xlabel("Time (s)")

    if units:
        ax.set_ylabel(f"{label} ({units})")
    else:
        ax.set_ylabel(label)

    if lower_constraint:
        ax.axhline(y=lower_constraint, color='r', linestyle='--', label=lower_constraint_label)
    if upper_constraint:
        ax.axhline(y=upper_constraint, color='g', linestyle='--', label=upper_constraint_label)
        
    ax.set_title(f"{label} vs. Time")
    ax.grid(True)
    ax.legend()
    # plt.show()


def plot_vec_madnitudes(x, y, z, label, units=None, lower_constraint=None, upper_constraint=None, lower_constraint_label=None, upper_constraint_label=None, sigma=15, ax=None):
    # create new figure
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
    vec = np.sqrt(x**2 + y**2 + z**2)
    time = np.linspace(0, sigma, len(vec))

    ax.plot(time, vec, label=label, color='blue', linewidth=1.3, marker='o', markersize=4)  # plot the data
    ax.set_xlabel("Time (s)")

    if units:
        ax.set_ylabel(f"{label} ({units})")
    else:
        ax.set_ylabel(label)

    if lower_constraint:
        ax.axhline(y=lower_constraint, color='r', linestyle='--', label=lower_constraint_label)
    if upper_constraint:
        ax.axhline(y=upper_constraint, color='g', linestyle='--', label=upper_constraint_label)

    ax.set_title(f"{label} vs. Time")
    ax.grid(True)
    ax.legend()
    # plt.show()




if __name__ == "__main__":

    # Example usage:
    file_path = "test/unit/trajectory_data.h5" # this is where the file will write
    # file_path = "jake_results/0/trajectory_data.h5" 
    # file_path = "jake_results/1/trajectory_data.h5" #i like
    # file_path = "jake_results/2/trajectory_data.h5" 
    # file_path = "jake_results/3/trajectory_data.h5" # i like
    # file_path = "jake_results/4/trajectory_data.h5"
    # file_path = "jake_results/5/trajectory_data.h5"
    # file_path = "jake_results/6/trajectory_data.h5"


    thrust_scale = 0.00004
    attitude_scale = 25
    # Read the data
    X_all, U_all, sigma_all = read_trajectory_data(file_path)
    sig = sigma_all[-1]

    # Plot the first trajectory
    # X_1 = X_all[0]
    # U_1 = U_all[0]
    # plot_trajectory(X_all[0], U_all[0], 1,  thrust_scale, attitude_scale)

    # plot the last trajectory
    X_last = X_all[-1]
    U_last = U_all[-1]
    plot_trajectory(X_last, U_last, len(X_all), thrust_scale, attitude_scale)

    # # plot all trajectories
    # for i in range(len(X_all)):
    #     plot_trajectory(X_all[i], U_all[i], i+1, thrust_scale, attitude_scale)


    iter_index = len(X_all) - 1

    # q_norm = np.sqrt(
    # X_all[iter_index][8, :]**2 +
    # X_all[iter_index][9, :]**2 +
    # X_all[iter_index][10, :]**2 +
    # X_all[iter_index][11, :]**2
    # )
    # # print(q_norm)
    # for i, q in enumerate(q_norm):
    #     X_all[iter_index][8:12, i] /= q
    # # X_all[iter_index][8:11, :] /= q_norm

    # # X_all[iter_index][8:11, :] /= q_norm
    # q_norm = np.sqrt(
    # X_all[iter_index][8, :]**2 +
    # X_all[iter_index][9, :]**2 +
    # X_all[iter_index][10, :]**2 +
    # X_all[iter_index][7, :]**2
    # )
    # # print(q_norm)
    # # # iter_index = 1

    # arg = 1 - 2*(X_all[iter_index][9, :]**2 - X_all[iter_index][10, :]**2)
    # for i, arg in enumerate(arg):
    #     if arg < -1:
    #         print(arg)
    #     elif arg > 1:
    #         print(arg)
    #     print(i, arg)

    # values = np.array([-0.0018703979885228293, -0.008111936695115203, -0.012287922923741354, 
    #                -0.014398470748599278, -0.015012139217765421, -0.014849214993926859, 
    #                -0.014501929047820852, -0.01363670695637376, -0.010975715571155197, 
    #                -0.0038132062843877993])
    # print("Input to arccos:", values)
    # print("Are values within [-1, 1]?", np.all((values >= -1) & (values <= 1)))
    # results = np.arccos(values)
    # print("arccos results:", results)
    # ## MASS 
    # # wet_mass = 30000.0
    # # dry_mass = 22000.0
    # # plot_state_var(X_all[iter_index][0, :], "Mass", "kg", lower_constraint=dry_mass, upper_constraint=wet_mass, lower_constraint_label="Dry Mass", upper_constraint_label="Wet Mass", sigma=sig)

    # # # X POS
    # # plot_state_var(X_all[iter_index][1, :], "X Position", "m", sigma=sig)

    # # # Y POS
    # # plot_state_var(X_all[iter_index][2, :], "Y Position", "m", sigma=sig)

    # # # Z POS
    # # plot_state_var(X_all[iter_index][3, :], "Z Position", "m", sigma=sig)

    # # X VEL
    # plot_state_var(X_all[iter_index][4, :], "X Velocity", "m/s", sigma=sig)

    # # Y VEL
    # plot_state_var(X_all[iter_index][5, :], "Y Velocity", "m/s", sigma=sig)

    # # Z VEL
    # plot_state_var(X_all[iter_index][6, :], "Z Velocity", "m/s", sigma=sig)

    # # Velocity magnitude
    # # plot_vec_madnitudes(X_all[iter_index][4, :], X_all[iter_index][5, :], X_all[iter_index][6, :], "Velocity", "m/s", sigma=sig)

    # ## QUATERNION q0
    # # plot_state_var(X_all[iter_index][7, :], "Quaternion q0", sigma=sig)

    # # # QUATERNION q1
    # # plot_state_var(X_all[iter_index][8, :], "Quaternion q1", sigma=sig)

    # # # QUATERNION q2
    # # plot_state_var(X_all[iter_index][9, :], "Quaternion q2", sigma=sig)

    # # # QUATERNION q3
    # # plot_state_var(X_all[iter_index][10, :], "Quaternion q3", sigma=sig)

    # # Conversion factor from radians to degrees
    RAD_TO_DEG = 180 / np.pi
    # # Convert angular velocities to degrees per second
    # # plot_state_var(X_all[iter_index][11, :] * RAD_TO_DEG, "Angular Velocity x", "deg/s", sigma=sig)
    # # plot_state_var(X_all[iter_index][12, :] * RAD_TO_DEG, "Angular Velocity y", "deg/s", sigma=sig)
    # # plot_state_var(X_all[iter_index][13, :] * RAD_TO_DEG, "Angular Velocity z", "deg/s", sigma=sig)
    
    # # plot_vec_madnitudes(X_all[iter_index][4, :], X_all[iter_index][5, :], X_all[iter_index][6, :], "Velocity", "m/s", sigma=sig)
    max_angular_velocity = 90 # deg/s
    # plot_vec_madnitudes(X_all[iter_index][11, :]* RAD_TO_DEG, X_all[iter_index][12, :]* RAD_TO_DEG, X_all[iter_index][13, :]* RAD_TO_DEG, "Angular Velocity", "deg/s", upper_constraint=max_angular_velocity, upper_constraint_label="Maximum Angular Velocity",sigma=sig)
    
    # ## THRUST 
    T_max = 800000.0
    T_min = 320000.0
    # plot_vec_madnitudes(U_all[iter_index][0, :], U_all[iter_index][1, :], U_all[iter_index][2, :], "Thrust", "N", lower_constraint=T_min, upper_constraint=T_max, lower_constraint_label="Min Thrust", upper_constraint_label="Max Thrust", sigma=sig)

    # # Tilt angle
    theta_max = 70.0
    cos_theta_max = 0.34202

    # theta_max = np.sqrt((1 - cos_theta_max) / 2)
    # this is index 10 and 11 in julia
    arg1 = 1 - 2*(X_all[iter_index][7, :]**2 + X_all[iter_index][8, :]**2) # q0 and q1
    arg2 = 1 - 2*(X_all[iter_index][8, :]**2 + X_all[iter_index][9, :]**2) # q1 and q2
    arg3 = 1 - 2*(X_all[iter_index][9, :]**2 + X_all[iter_index][10, :]**2) # q2 and q3
    arg4 = 1 - 2*(X_all[iter_index][10, :]**2 + X_all[iter_index][7, :]**2) # q3 and q0
    tilt1 = np.arccos(arg1) * RAD_TO_DEG
    tilt2 = np.arccos(arg2) * RAD_TO_DEG
    tilt3 = np.arccos(arg3) * RAD_TO_DEG
    tilt4 = np.arccos(arg4) * RAD_TO_DEG

    # make a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    plot_state_var(tilt1, "Tilt Angle 1", "deg", upper_constraint=theta_max, upper_constraint_label="Max Tilt Angle", sigma=sig, ax=axs[0, 0])
    plot_state_var(tilt2, "Tilt Angle 2", "deg", upper_constraint=theta_max, upper_constraint_label="Max Tilt Angle", sigma=sig, ax=axs[0, 1])
    plot_state_var(tilt3, "Tilt Angle 3", "deg", upper_constraint=theta_max, upper_constraint_label="Max Tilt Angle", sigma=sig, ax=axs[1, 0])
    plot_state_var(tilt4, "Tilt Angle 4", "deg", upper_constraint=theta_max, upper_constraint_label="Max Tilt Angle", sigma=sig, ax=axs[1, 1])



    arg = np.clip(arg1, -1, 1)
    tilt_angle = np.arccos(arg)
    tilt_angle = compute_tilt_angles(X_all[iter_index][7:11, :])
    # tilt_angle = X_all[iter_index][9, :]**2 + X_all[iter_index][10, :]**2
    tilt_angle = tilt_angle * RAD_TO_DEG

    # plot_state_var(tilt_angle, "Tilt Angle", "deg", upper_constraint=theta_max, upper_constraint_label="Max Tilt Angle", sigma=sig)

    # ## Plot the up position on the y axis and the east position on the x axis
    # plot_trajectory_2D(X_all, U_all, iter_index, thrust_scale, attitude_scale, (1, 3), (0, 2))
    # plot_trajectory_2D(X_all, U_all, iter_index, thrust_scale, attitude_scale, (2, 3), (1, 2))
    # plot_trajectory_2D(X_all, U_all, iter_index, thrust_scale, attitude_scale, (1, 2), (0, 1))

    # # ## Plot the north position on the y axis and the east position on the x axis
    # # fig = plt.figure(figsize=(8, 6))
    # # ax = fig.add_subplot(111)
    # # ax.plot(X_all[iter_index][1, :], X_all[iter_index][2, :], label="Trajectory", color='blue', linewidth=1.3, marker='o')  # plot the data
    # # ax.set_xlabel("East Position (m)")
    # # ax.set_ylabel("North Position (m)")
    # # ax.set_title("North vs East Position")
    # # ax.grid(True)
    # # ax.legend()
    
    # # # ## Plot the up position on the y axis and the north position on the x axis
    # # fig = plt.figure(figsize=(8, 6))
    # # ax = fig.add_subplot(111)
    # # ax.plot(X_all[iter_index][2, :], X_all[iter_index][3, :], label="Trajectory", color='blue', linewidth=1.3, marker='o')  # plot the data
    # # ax.set_xlabel("North Position (m)")
    # # ax.set_ylabel("Up Position (m)")
    # # ax.set_title("Up vs North Position")
    # # ax.grid(True)
    # # ax.legend()
    
    # # ## Plot the three components of the angular velocity vs time
    # fig = plt.figure(figsize=(8, 6))
    # ax = fig.add_subplot(111)
    # time = np.linspace(0, sig, len(X_all[iter_index][11, :]))
    # ax.plot(time, X_all[iter_index][11, :] * RAD_TO_DEG, label="Angular Velocity x", color='blue', linewidth=1.3, marker='o')  # plot the data
    # ax.plot(time, X_all[iter_index][12, :] * RAD_TO_DEG, label="Angular Velocity y", color='red', linewidth=1.3, marker='o')  # plot the data
    # ax.plot(time, X_all[iter_index][13, :] * RAD_TO_DEG, label="Angular Velocity z", color='green', linewidth=1.3, marker='o')  # plot the data
    # ax.set_xlabel("Time step")
    # ax.set_ylabel("Angular Velocity (deg/s)")
    # ax.set_title("Angular Velocity vs Time")
    # ax.grid(True)
    # ax.legend()


    # define a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    # ax_3d = fig.add_subplot(2, 2, 4, projection='3d')

    plot_trajectory_2D(X_all, U_all, iter_index, thrust_scale, attitude_scale, (1, 3), (0, 2), ax=axs[0, 0])
    plot_trajectory_2D(X_all, U_all, iter_index, thrust_scale, attitude_scale, (2, 3), (1, 2), ax=axs[0, 1])
    plot_trajectory_2D(X_all, U_all, iter_index, thrust_scale, attitude_scale, (1, 2), (0, 1), ax=axs[1, 0])
    # plot_trajectory(X_all[iter_index], U_all[iter_index], iter_index, thrust_scale, attitude_scale, ax=ax_3d)
    plot_vec_madnitudes(X_all[iter_index][4, :], X_all[iter_index][5, :], X_all[iter_index][6, :], "Velocity", "m/s", sigma=sig, ax=axs[1, 1])
    # plot_vec_madnitudes(U_all[iter_index][0, :], U_all[iter_index][1, :], U_all[iter_index][2, :], "Thrust", "N", lower_constraint=T_min, upper_constraint=T_max, lower_constraint_label="Min Thrust", upper_constraint_label="Max Thrust", sigma=sig, ax=axs[1, 1])
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    # plot the thrust and the tilt angle
    plot_vec_madnitudes(U_all[iter_index][0, :], U_all[iter_index][1, :], U_all[iter_index][2, :], "Thrust", "N", lower_constraint=T_min, upper_constraint=T_max, lower_constraint_label="Min Thrust", upper_constraint_label="Max Thrust", sigma=sig, ax=axs[0])
    plot_state_var(tilt_angle, "Tilt Angle", "deg", upper_constraint=theta_max, upper_constraint_label="Max Tilt Angle", sigma=sig, ax=axs[1])
    
    plt.tight_layout()
    plt.show()


