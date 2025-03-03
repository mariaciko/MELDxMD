#! /usr/bin/env python

#ATOM      3  CA  ARG     1       3.970   2.846  -0.000  1.00  0.00
#ATOM    142 HH22 ARG     6      26.753   9.204  -4.004  1.00  0.00

import mdtraj as md
import numpy as np

traj = md.load("renumbered.pdb")
protein_chain = traj.topology.select('chainid == 0')
peptide_chain = traj.topology.select('chainid == 1')

initial_com = md.compute_center_of_mass(traj.atom_slice(protein_chain))
translation_vector = -initial_com  # Negate to shift to origin
peptide_translation_vector = np.array([0, 0, 4]) # shift peptide 4 nm away from protein
traj.xyz[:, protein_chain, :] += translation_vector # the : means all frames, all atoms in the protein chain, and : means all coordinates (x, y, z)
traj.xyz[:, peptide_chain, :] += peptide_translation_vector

# traj.save('pep_shifted.pdb')
traj.save("TEMPLATES/pep_shifted.pdb")
final_com = md.compute_center_of_mass(traj.atom_slice(protein_chain))

print(f'Initial COM protein: {initial_com}')
print(f'Final COM protein: {final_com}')
