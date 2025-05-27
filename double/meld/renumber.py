import mdtraj as md

# Load trajectory and topology
traj = md.load("stacked.pdb")

# Identify protein and peptide residues
protein_residues = [res for res in traj.topology.residues if res.index <= 84]
peptide1_residues = [res for res in traj.topology.residues if res.index > 84]  
peptide2_residues = [res for res in traj.topology.residues if res.index > 96]

# Renumber residues: Protein 0-84, Peptide1 85-96, Peptide2 97-108
for i, res in enumerate(protein_residues):
    res.resSeq = i + 0 

for i, res in enumerate(peptide1_residues):
    res.resSeq = i + 85 

for i, res in enumerate(peptide2_residues):
    res.resSeq = i + 97

# Save the modified topology
traj.save_pdb("renumbered.pdb")

print("Residue renumbering completed. Protein: 0-84, Peptide1: 85-96, Peptide2: 97-108.")
