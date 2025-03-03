import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory and topology files
traj0 = md.load('traj.dcd', top='topology.pdb')
traj29 = md.load("traj29.dcd",top="topology29.pdb")

def compute_com_distance(traj_file: str, top_file: str, label: str):

    traj = md.load(traj_file, top=top_file)

    protein_indices = traj.top.select('resid 0 to 84')
    peptide_indices = traj.top.select('resid 85 to 96')

    protein_traj = traj.atom_slice(protein_indices)
    peptide_traj = traj.atom_slice(peptide_indices)

    com_protein = md.compute_center_of_mass(protein_traj)
    com_peptide = md.compute_center_of_mass(peptide_traj)

    # Calculate the distance between the centers of mass
    distances = np.linalg.norm(com_protein - com_peptide, axis=1) # linalg.norm computes the Euclidean norm (magnitude); axis=1 means we are calculating the norm along the second axis because the first axis is the frame number

    # Create a time array for the x-axis
    num_steps = traj.n_frames # Number of frames (steps)
    timestep_fs = 4.5  # Timestep in femtoseconds (4.5 fs = 4.5e-3 ps)
    total_time_ps = num_steps * timestep_fs / 1000  # Convert to picoseconds
    time = np.linspace(0, total_time_ps, num_steps)  # Create a time array from 0 to total time

    return time, distances, label

# File paths for replica 0 and replica 29
traj_0, top_0 = 'traj.dcd', 'topology.pdb'
traj_29, top_29 = 'traj29.dcd', 'topology29.pdb'

# Compute distances for both replicas
time_0, distances_0, label_0 = compute_com_distance(traj_0, top_0, 'R0: Distance between COMs')
time_29, distances_29, label_29 = compute_com_distance(traj_29, top_29, 'R29: Distance between COMs')

# Plot distance as a function of time for both replicas
plt.figure(figsize=(10, 6))
plt.plot(time_0, distances_0, label=label_0, color='blue')
plt.plot(time_29, distances_29, label=label_29, color='red')
plt.xlim(-5, max(time_0[-1], time_29[-1]) +5) # the max function is used to ensure that the x-axis limit is set to the maximum value between time_0 and time_29
plt.xticks(np.arange(0, max(time_0[-1], time_29[-1]) +5, 10))  # Set x-ticks at 0, 30, 60, 90 ps
plt.xlabel('Time (ps)')
plt.ylabel('Distance (nm)')
plt.title('Distance between Center of Mass of Protein and Peptide over Time')
plt.grid()
plt.legend()
plt.tight_layout()

# Save plot
plt.savefig('COM_distance.png', dpi=300)
plt.close()
