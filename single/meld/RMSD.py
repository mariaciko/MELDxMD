import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

# Load trajectory and topology files
pdb = md.load("renumbered.pdb")
pdb1 = pdb.top.select("name CA and resid 0 to 84")
pdb2 = pdb.top.select("name CA and resid 85 to 96")
traj_data = {"Replica 0": ("traj.dcd", "topology.pdb"),
             "Replica 29": ("traj29.dcd", "topology29.pdb")}

rmsd_results = {}
time_results = {}

for label, (traj_file, top_file) in traj_data.items():
    traj = md.load(traj_file, top=top_file)
    
    # Define atom selections
    li = traj.top.select('name CA and resid 85 to 96')
    re = traj.top.select("name CA and resid 0 to 84")

    rmsd_re = md.rmsd(traj, pdb, atom_indices=re, ref_atom_indices=pdb1)
    rmsd_li = md.rmsd(traj, pdb, atom_indices=li, ref_atom_indices=pdb2)
    
    #LRMSD
    traj.superpose(pdb, atom_indices=re, ref_atom_indices=pdb1)
    reference_frame = pdb.xyz[0, pdb2]
    # print(len(reference_frame1))
    lrmsd = []

    for frame_index in range(0, len(traj)):
    # Extract the coordinates of the current frame for the specified atoms
        current_coordinates = traj.xyz[frame_index, li]

        # Calculate the squared differences between the current frame and the reference frame
        squared_diff = ((current_coordinates - reference_frame) ** 2).sum()

        # Calculate the RMSD using the formula (square root of the mean squared difference)
        rmsd = (squared_diff / len(current_coordinates)) ** 0.5

        # Append the RMSD value to the list
        lrmsd.append(rmsd)

    # Store results
    rmsd_results[label] = {"Ligand 1": rmsd_li, "Receptor": rmsd_re, "LRMSD ": lrmsd}

    # Define time axis
    num_frames = traj.n_frames
    timestep_fs = 4.5 
    total_time_ps = num_frames * timestep_fs / 1000  # Convert to ps
    time_results[label] = np.linspace(0, total_time_ps, num_frames)


# Create subplots (2 rows, 1 column)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True, sharey=True)

replica_labels = list(traj_data.keys())

# Plot RMSD for each replica
for ax, (label, rmsds) in zip(axes, rmsd_results.items()):
    mean_rmsd_li = np.mean(rmsds["Ligand 1"])
    mean_rmsd_re = np.mean(rmsds["Receptor"])
    mean_rmsd_lr = np.mean(rmsds["LRMSD "])

    ax.plot(time_results[label], rmsds["Ligand 1"], label="Ligand", color='r')
    ax.plot(time_results[label], rmsds["Receptor"], label="Receptor", color='b')
    ax.plot(time_results[label], rmsds["LRMSD "], label="LRMSD", color='y')

    ax.text(time_results[label][-1]+14, mean_rmsd_li, f"Mean:{mean_rmsd_li:.2f}", color='r')
    ax.text(time_results[label][-1]+14, mean_rmsd_re, f"Mean:{mean_rmsd_re:.2f}", color='b')
    ax.text(time_results[label][-1]+14, mean_rmsd_lr, f"Mean:{mean_rmsd_lr:.2f}", color='y')

    ax.axhline(y=mean_rmsd_li, color='black', linestyle='--', label="Ligand Mean")
    ax.axhline(y=mean_rmsd_re, color='black', linestyle='--', label="Receptor Mean")
    ax.axhline(y=mean_rmsd_lr, color='black', linestyle='--', label="LRMSD Mean")

    ax.set_ylabel("RMSD (nm)")
    ax.set_title(f"{label}: MDM2-P53 RMSD over Time")
    ax.grid()
    ax.legend(["RMSD-L", "RMSD-R", "LRMSD"], loc='upper right')

# X-axis label for the bottom subplots
for ax in axes[1:]:
    ax.set_xlabel("Time (ps)")

# plt.tight_layout()
plt.savefig("RMSD_subplots.png", dpi=300)
plt.close()