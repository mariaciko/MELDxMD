import numpy as np
import mdtraj as md
import os
import matplotlib.pyplot as plt

# Load topology and check atom count
topology = md.load_pdb("topology.pdb")
print("Number of atoms in topology:", topology.n_atoms)
trajectory = md.load("traj.dcd", top="topology.pdb")
print("Number of atoms in trajectory:", trajectory.n_atoms)
print("Number of frames in trajectory:", trajectory.n_frames)

# Define protein and peptide residue indices
protein_indices = trajectory.topology.select('resid 0 to 84')
peptide_indices = trajectory.topology.select('resid 85 to 96')

# Directory to save plots
output_dir = "COM_distance_plots"
os.makedirs(output_dir, exist_ok=True)

# Compute COM distance for each runner
num_runners = 30
for runner in range(num_runners):
    runner_frames = trajectory[runner::num_runners]
    
    # Compute COM for protein and peptide for each frame
    protein_com = md.compute_center_of_mass(runner_frames.atom_slice(protein_indices))
    peptide_com = md.compute_center_of_mass(runner_frames.atom_slice(peptide_indices))
    
    # Calculate the COM distance
    com_distances = np.linalg.norm(protein_com - peptide_com, axis=1)
    
    # Smoothing using a moving average
    window_size = 10  # You can adjust this
    smoothed_distances = np.convolve(com_distances, np.ones(window_size)/window_size, mode='valid')
    
    # Calculate mean and standard deviation of the COM distances
    mean_distance = np.mean(smoothed_distances)
    std_distance = np.std(smoothed_distances)

    # Set thresholds based on mean and standard deviation
    docked_threshold = mean_distance - std_distance  # Example: one standard deviation below the mean
    unbound_threshold = mean_distance + std_distance  # Example: one standard deviation above the mean

    # Create a mask for docked and unbound states
    state_colors = ['red' if distance > unbound_threshold else 'green' for distance in smoothed_distances]
    
    # Define total number of MELD steps
    total_steps = 7848
    timestep = 4.5  # in fs
    num_frames = len(smoothed_distances)  # Number of frames after smoothing

    # Create x-axis values corresponding to MELD steps
    meld_steps = np.arange(num_frames) * (total_steps / (num_frames + window_size - 1))

    # Plot the COM distance for each runner
    plt.figure(figsize=(10, 6))
    plt.plot(meld_steps, smoothed_distances, label=f'Runner {runner + 1}', color='blue')
    
    # Highlight docked (green) and unbound (red) states
    for i in range(len(smoothed_distances) - 1):
        plt.plot([meld_steps[i], meld_steps[i + 1]], smoothed_distances[i:i+2], color=state_colors[i])
    
    plt.axhline(y=docked_threshold, color='green', linestyle='--', label='Docked Threshold')
    plt.axhline(y=unbound_threshold, color='red', linestyle='--', label='Unbound Threshold')
    
    plt.xlabel("MELD Steps")
    plt.ylabel("COM Distance (nm)")
    plt.title(f"COM Distance - Runner {runner + 1} (Start of sim)")
    plt.legend()
    
    # Save the plot
    plt.savefig(f"{output_dir}/COM_distance_runner_{runner + 1}.png")
    plt.close()

print(f"Plots saved in '{output_dir}' directory.")
