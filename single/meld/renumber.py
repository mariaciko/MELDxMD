import mdtraj as md

# Load trajectory and topology
traj = md.load("../native.pdb")

# Identify protein and peptide residues
protein_residues = [res for res in traj.topology.residues if res.index < 85]
peptide_residues = [res for res in traj.topology.residues if res.index >= 85]  

# Renumber residues: Protein from 1 to 85, Peptide from 86 to 97
for i, res in enumerate(protein_residues):
    res.resSeq = i + 0  # Start from 0

for i, res in enumerate(peptide_residues):
    res.resSeq = i + 85  # Start from 85

# Save the modified topology
traj.save_pdb("renumbered.pdb")

print("Residue renumbering completed. Protein: 0-84, Peptide: 85-96.")
