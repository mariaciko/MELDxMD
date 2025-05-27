import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

# Load trajectory and topology files
pdb = md.load("renumbered.pdb")
pdb1 = pdb.top.select("name CA and resid 0 to 84")
pdb2 = pdb.top.select("name CA and resid 85 to 96")
pdb3 = pdb.top.select("name CA and resid 97 to 108")
traj_data = {"Replica 0": ("traj.dcd", "topology.pdb"),
             "Replica 29": ("traj29.dcd", "topology29.pdb")}

rmsd_results = {}
time_results = {}

for label, (traj_file, top_file) in traj_data.items():
    traj = md.load(traj_file, top=top_file)
    
    # Define atom selections
    li1 = traj.top.select('name CA and resid 85 to 96')
    li2 = traj.top.select('name CA and resid 97 to 108')
    re = traj.top.select("name CA and resid 0 to 84")

    rmsd_re1 = md.rmsd(traj, pdb, atom_indices=re, ref_atom_indices=pdb1)
    rmsd_li1 = md.rmsd(traj, pdb, atom_indices=li1, ref_atom_indices=pdb2)
    rmsd_li2 = md.rmsd(traj, pdb, atom_indices=li2, ref_atom_indices=pdb3)
    
    #LRMSD1
    traj.superpose(pdb, atom_indices=re, ref_atom_indices=pdb1)
    reference_frame1 = pdb.xyz[0, pdb2]
    # print(len(reference_frame1))
    lrmsd1 = []

    for frame_index in range(0, len(traj)):
    # Extract the coordinates of the current frame for the specified atoms
        current_coordinates = traj.xyz[frame_index, li1]

        # Calculate the squared differences between the current frame and the reference frame
        squared_diff = ((current_coordinates - reference_frame1) ** 2).sum()

        # Calculate the RMSD using the formula (square root of the mean squared difference)
        rmsd = (squared_diff / len(current_coordinates)) ** 0.5

        # Append the RMSD value to the list
        lrmsd1.append(rmsd)
    
    
    #LRMSD2
    traj.superpose(pdb, atom_indices=re, ref_atom_indices=pdb1)
    reference_frame2 = pdb.xyz[0, pdb2]
    print(len(reference_frame2))
    lrmsd2 = []

    for frame_index in range(0, len(traj)):
    # Extract the coordinates of the current frame for the specified atoms
        current_coordinates = traj.xyz[frame_index, li2]

        # Calculate the squared differences between the current frame and the reference frame
        squared_diff = ((current_coordinates - reference_frame2) ** 2).sum(axis=1)

        # Calculate the RMSD using the formula (square root of the mean squared difference)
        rmsd = (squared_diff / len(current_coordinates)) ** 0.5

        # Append the RMSD value to the list
        lrmsd2.append(rmsd)

    # Store results
    rmsd_results[label] = {"Ligand 1": rmsd_li1, "Receptor": rmsd_re1, "LRMSD 1": lrmsd1,
                           "Ligand 2": rmsd_li2,"LRMSD 2": lrmsd2}

    # Define time axis
    num_frames = traj.n_frames
    timestep_fs = 4.5 
    total_time_ps = num_frames * timestep_fs / 1000  # Convert to ps
    time_results[label] = np.linspace(0, total_time_ps, num_frames)


# Create subplots (2 rows, 2 columns)
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8), sharex=True, sharey=True)

replica_labels = list(traj_data.keys())
ligand_labels = ["Ligand 1", "Ligand 2"]

# Plot RMSD for each replica-ligand pair
for i, rep_label in enumerate(replica_labels):
    for j, lig_label in enumerate(ligand_labels):
        ax = axes[i, j]
        ax.plot(time_results[rep_label], rmsd_results[rep_label][lig_label], label="RMSD-L", color='r')
        ax.plot(time_results[rep_label], rmsd_results[rep_label]["Receptor"], label="RMSD-R", color='b')
        ax.plot(time_results[rep_label], rmsd_results[rep_label][f"LRMSD {j+1}"], label=f"LRMSD {lig_label[-1]}", color='y')
        ax.set_ylabel("RMSD (nm)")
        ax.set_title(f"{rep_label}: {lig_label} RMSD over Time")
        ax.grid()
        ax.legend(["RMSD-L", "RMSD-R", f"LRMSD {lig_label[-1]}"], loc='upper right')

# X-axis label for the bottom subplots
for ax in axes[-1, :]:
    ax.set_xlabel("Time (ps)")

plt.tight_layout()
plt.savefig("RMSD_subplots.png", dpi=300)
