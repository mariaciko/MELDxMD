import numpy as np
import mdtraj as md

# Load PDB structure
pdb = md.load('renumbered.pdb')

# Select CA atoms of the protein
CA_protein_indices = pdb.top.select("name == CA and chainid == 0")
print("Protein CA indices: ", CA_protein_indices)
CA_protein = pdb.atom_slice(CA_protein_indices)
print("Number of protein CA atoms: ", CA_protein.n_atoms)

# Select CA atoms of peptide
CA_peptide_indices = pdb.top.select("name == CA and chainid == 1")
print("Peptide CA indices: ", CA_peptide_indices)
CA_peptide = pdb.atom_slice(CA_peptide_indices)
print("Peptide CA atoms: ", CA_peptide.n_atoms)


# Select CA pairs
pairs = pdb.top.select_pairs(CA_protein_indices, CA_peptide_indices)
print("CA Pairs: ", pairs)

# Compute distances between the pairs
dist = md.compute_distances(pdb, pairs)[0] # returns a tuple, we only need the first element
print("Distance bw protein & peptide1 CA atoms: ", dist)

# Apply nonbonded cutoff distance
cutoff_distance = 0.9
filtered_indices = (dist < cutoff_distance)
print("Filtered indices: ", filtered_indices)
print("Shape of filtered indices: ", np.shape(filtered_indices))
print("Number of filtered indices: ", np.sum(filtered_indices))

# Filter the distances
filtered_dist = dist[filtered_indices]
print("Filtered distances: ", filtered_dist)
print("Shape: ", np.shape(filtered_dist))
print("Number: ", np.sum(filtered_dist))

# Filter the pairs
contact_pairs = pairs[filtered_indices]
print("Peptide contact pairs: ", contact_pairs)


for i, contact in enumerate(contact_pairs):
    print(pdb.top.atom(contact[0]).residue.resSeq, pdb.top.atom(contact[0]).residue.name, "CA", 
          pdb.top.atom(contact[1]).residue.resSeq, pdb.top.atom(contact[1]).residue.name, "CA",
          filtered_dist[i])



with open("protein_pep_all.dat", "w") as fout:
    for i, contact in enumerate(contact_pairs):
        r1 = pdb.top.atom(contact[0]).residue.resSeq # residue number of the first atom in the contact pair
        r2 = pdb.top.atom(contact[1]).residue.resSeq
        n1 = pdb.top.atom(contact[0]).residue.name
        n2 = pdb.top.atom(contact[1]).residue.name
        string = "{:d} {:s} {:} {:d} {:s} {:} {:3.5f}\n".format(r1, n1, "CA", r2, n2, "CA", filtered_dist[i])
        fout.write(string)
        fout.write("\n")


with open ("pseudobonds_complex.pb", "w") as fout:
    for i, contact in enumerate(contact_pairs):
        string = "/A:{:d}@CA /B:{:d}@CA\n".format(pdb.top.atom(contact[0]).residue.resSeq,
                                                  pdb.top.atom(contact[1]).residue.resSeq)
        fout.write(string)