import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

# # Load the trajectory and topology files
# topology0 = md.load_pdb("topology.pdb")
# traj0 = md.load("traj.dcd", top="topology.pdb")
# topology29 = md.load_pdb("topology29.pdb")
# traj29 = md.load("traj29.dcd", top="topology29.pdb")

# for traj_file in [traj0, traj29]:
#     print("Number of atoms:", traj_file.n_atoms)
#     li = traj_file.top.select('resid 85 to 97')
#     print("Number of ligand atoms:", len(li))

#     traj_file.superpose(traj_file[0], atom_indices=li)
#     rmsd_li = md.rmsd(traj_file,traj_file[0],atom_indices=li)
#     # print("RMSD of ligand:", rmsd_li)

#     reference_frame = traj_file.xyz[0, li]
#     rmsd_values = []



# for frame_index in range(1, len(traj)):
#     # Extract the coordinates of the current frame for the specified atoms
#     current_coordinates = traj.xyz[frame_index, li]

#     # Calculate the squared differences between the current frame and the reference frame
#     squared_diff = ((current_coordinates - reference_frame) ** 2).sum()

#     # Calculate the RMSD using the formula (square root of the mean squared difference)
#     rmsd = (squared_diff / len(current_coordinates)) ** 0.5

#     # Append the RMSD value to the list
#     rmsd_values.append(rmsd)

# # Time in ps for each sampled frame
# timestep_fs = 4.5
# timestep_ps = timestep_fs / 1000  # Convert to ps
# time_ps = np.arange(0, len(rmsd_li) * timestep_ps, 1000 * timestep_ps)

# # Ensure time ends at 90 ps
# print(f"Total simulation time: {time_ps[-1]} ps")

# plt.plot(rmsd_li, 'b', label='Ligand')


# # Save every 1000th frame as a PDB file for visualization
# for i in range(0, len(traj), 1000):
#     frame = traj[i]  # Select the i-th frame (every 1000th frame)
#     filename = f"0frame_{i}.pdb"  # PDB filename for each frame
#     frame.save_pdb(filename)  # Save the frame as a PDB file
#     print(f"Saved {filename}")



# Load trajectory and topology files
traj_data = {
    "Replica 0": ("traj.dcd", "topology.pdb"),
    "Replica 29": ("traj29.dcd", "topology29.pdb"),
}

rmsd_results = {}
time_results = {}

for label, (traj_file, top_file) in traj_data.items():
    traj = md.load(traj_file, top=top_file)
    
    # Select CA ligand atoms
    ligand_indices = traj.top.select('resid 85 to 96')
    CA_ligand_indices = traj.top.select('resid 85 to 96 and name CA')
    print(f"{label} - Number of CA ligand atoms: {len(CA_ligand_indices)}")

    # Align trajectory to the first frame using ligand atoms
    traj.superpose(traj[0], atom_indices=ligand_indices)

    # Compute RMSD
    rmsd_ligand = md.rmsd(traj, traj[0], atom_indices=ligand_indices)
    
    # Store results
    rmsd_results[label] = rmsd_ligand

    # Define time axis
    num_frames = traj.n_frames
    timestep_fs = 4.5 
    total_time_ps = num_frames * timestep_fs / 1000  # Convert to ps
    time_results[label] = np.linspace(0, total_time_ps, num_frames)

# Create subplots (2 rows, 1 column)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)

# Plot RMSD for each replica in separate subplots
for ax, (label, rmsd_ligand) in zip(axes, rmsd_results.items()):
    ax.plot(time_results[label], rmsd_ligand, label=label, color='b')
    ax.set_ylabel("Ligand RMSD (nm)")
    ax.set_title(f"{label}: MDM2-P53 Ligand RMSD over Time")
    ax.grid()
    ax.legend()

# X-axis label for the last subplot
axes[-1].set_xlabel("Time (ps)")

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig("RMSD_subplots.png", dpi=300)
plt.close()
