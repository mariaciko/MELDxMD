import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory and topology
traj = md.load("traj29.dcd", top="topology29.pdb")

# Select atoms for alignment and peptide COM calculations
protein_atoms = traj.topology.select('resid 0 to 84')
peptide_atoms = traj.topology.select('resid 85 to 96')

# Align protein to the first frame
traj.superpose(traj, atom_indices=protein_atoms)

# Initialize list to store COM of peptide in each frame
peptide_com_list = []

# Calculate COM of peptide for each frame
for i, frame in enumerate(traj):
    if i % 4 == 0:  # Save COM every 4 frames
        peptide_com = md.compute_center_of_mass(frame.atom_slice(peptide_atoms))
        peptide_com = np.squeeze(peptide_com)  # Remove the extra dimension if present
        if peptide_com.shape == (3,):  # Ensure the output is a 3D coordinate
            peptide_com_list.append(peptide_com)
        else:
            print(f"Unexpected shape in COM calculation after squeeze: {peptide_com.shape}")

# Convert to numpy array for easier manipulation
peptide_com_array = np.array(peptide_com_list)

# Save COM data to an XYZ file for ChimeraX
with open('peptide_com.xyz', 'w') as file:
    file.write(f"{len(peptide_com_array)}\n")
    file.write("COM coordinates for every 4th frame\n")
    for com in peptide_com_array:
        file.write(f"H {com[0]} {com[1]} {com[2]}\n")

# Convert COM coordinates to spherical coordinates relative to protein center
protein_com = np.mean(traj.xyz[:, protein_atoms, :], axis=1)
relative_com = peptide_com_array - protein_com[:, np.newaxis]

# Convert Cartesian to spherical coordinates
r = np.linalg.norm(relative_com, axis=1) # radial distance 
theta = np.arctan2(relative_com[:, 1], relative_com[:, 0])  # Azimuthal angle
phi = np.arccos(relative_com[:, 2] / r)  # Polar angle

# Plot and save histograms
plt.figure(figsize=(8, 5))
plt.hist(theta, bins=50, density=True)
plt.xlabel(r"$\theta$ (radians)")
plt.ylabel(r"$p(\theta)$")
plt.title(r"$p(\theta)$ vs $\theta$")
plt.savefig("p_theta_vs_theta.png")
plt.close()

plt.figure(figsize=(8, 5))
plt.hist(phi, bins=50, density=True)
plt.xlabel(r"$\phi$ (radians)")
plt.ylabel(r"$p(\phi)$")
plt.title(r"$p(\phi)$ vs $\phi$")
plt.savefig("p_phi_vs_phi.png")
plt.close()

plt.figure(figsize=(8, 5))
plt.hist(r, bins=50, density=True)
plt.xlabel(r"$r$ (nm)")
plt.ylabel(r"$p(r)$")
plt.title(r"$p(r)$ vs $r$")
plt.savefig("p_r_vs_r.png")
plt.close()

print("Plots saved as 'p_theta_vs_theta.png', 'p_phi_vs_phi.png', and 'p_r_vs_r.png'")