import mdtraj as md
import numpy as np

# Load trajectory and topology
traj = md.load('traj.dcd', top='topology.pdb')

# Select atoms (adjust for your system)
protein = traj.top.select("resid 0 to 84")
peptide1 = traj.top.select("resid 85 to 96")
peptide2 = traj.top.select("resid 97 to 108")

# Compute COM distances
com_distance_1 = md.compute_center_of_mass(traj.atom_slice(protein)) - md.compute_center_of_mass(traj.atom_slice(peptide1))
com_distance_2 = md.compute_center_of_mass(traj.atom_slice(protein)) - md.compute_center_of_mass(traj.atom_slice(peptide2))

R = 8.314  # J/(mol*K)
T = 300  # Kelvin

P1 = np.sum(com_distance_1 < 0.9) / len(com_distance_1) # Fraction of time peptide 1 is within 0.9 nm of the protein
P2 = np.sum(com_distance_2 < 0.9) / len(com_distance_2)
print(f"Peptide 1 binding probability: {P1:.2f}")
print(f"Peptide 2 binding probability: {P2:.2f}")

delta_G = -R * T * np.log(P2 / P1)  # in Joules
delta_G_kcal = delta_G / 4184  # Convert to kcal/mol


# Save results
with open('binding_free_energy.txt', 'w') as f:
    f.write(f"Relative binding free energy ΔG: {delta_G_kcal:.2f} kcal/mol\n")
