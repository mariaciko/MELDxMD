import numpy as np
import mdtraj as md

# Load PDB structure
pdb = md.load('renumbered.pdb')

# Select CA atoms of the protein
CA_protein_indices = pdb.top.select("name == CA and chainid == 0")
print("Protein CA indices: ", CA_protein_indices)
CA_protein = pdb.atom_slice(CA_protein_indices)
print("Number of protein CA atoms: ", CA_protein.n_atoms)

# Select CA atoms of peptides
CA_peptide_indices1 = pdb.top.select("name == CA and chainid == 1")
print("Peptide1 CA indices: ", CA_peptide_indices1)
CA_peptide1 = pdb.atom_slice(CA_peptide_indices1)
print("Peptide1 CA atoms: ", CA_peptide1.n_atoms)

CA_peptide_indices2 = pdb.top.select("name == CA and chainid == 2")
print("Peptide2 CA indices: ", CA_peptide_indices2)
CA_peptide2 = pdb.atom_slice(CA_peptide_indices2)
print("Peptide2 CA atoms: ", CA_peptide2.n_atoms)

# Select CA pairs
pairs1 = pdb.top.select_pairs(CA_protein_indices, CA_peptide_indices1)
print("CA Pairs1: ", pairs1)
pairs2 = pdb.top.select_pairs(CA_protein_indices, CA_peptide_indices2)
print("CA Pairs2: ", pairs2)

# Compute distances between the pairs
dist1 = md.compute_distances(pdb, pairs1)[0] # returns a tuple, we only need the first element
print("Distance bw protein & peptide1 CA atoms: ", dist1)
dist2 = md.compute_distances(pdb, pairs2)[0] # returns a tuple, we only need the first element
print("Distance bw protein & peptide2 CA atoms: ", dist2)

# Apply nonbonded cutoff distance
cutoff_distance = 0.9
filtered_indices1 = (dist1 < cutoff_distance)
print("Filtered indices: ", filtered_indices1)
print("Shape of filtered indices: ", np.shape(filtered_indices1))
print("Number of filtered indices: ", np.sum(filtered_indices1))

print("\n")

filtered_indices2 = (dist2 < cutoff_distance)
print("2nd: Filtered indices", filtered_indices2)
print("2nd: Shape of the filtered indices", np.shape(filtered_indices2))
print("2nd: Number of the filtered indices", np.sum(filtered_indices2))


# Filter the distances
filtered_dist1 = dist1[filtered_indices1]
print("Filtered distances: ", filtered_dist1)
print("Shape: ", np.shape(filtered_dist1))
print("Number: ", np.sum(filtered_dist1))

print("\n")

filtered_dist2 = dist2[filtered_indices2]
print("Second peptide: Filtered distances", filtered_dist2)
print("Second peptide: Shape of the filtered distances", np.shape(filtered_dist2))
print("Second peptide: Number of the filtered distances", np.sum(filtered_dist2))


# Filter the pairs
contact_pairs1 = pairs1[filtered_indices1]
print("Peptide contact pairs: ", contact_pairs1)
contact_pairs2 = pairs2[filtered_indices2]
print("Second peptide contact pairs", contact_pairs2)


for i, contact in enumerate(contact_pairs1):
    print(pdb.top.atom(contact[0]).residue.resSeq, pdb.top.atom(contact[0]).residue.name, "CA", 
          pdb.top.atom(contact[1]).residue.resSeq, pdb.top.atom(contact[1]).residue.name, "CA",
          filtered_dist1[i])

for i, contact in enumerate(contact_pairs2):
    print(pdb.top.atom(contact[0]).residue.resSeq, pdb.top.atom(contact[0]).residue.name, "CA", 
          pdb.top.atom(contact[1]).residue.resSeq, pdb.top.atom(contact[1]).residue.name, "CA",
          filtered_dist2[i])



with open("protein_pep_all.dat", "w") as fout:
    for i, contact in enumerate(contact_pairs1):
        r1 = pdb.top.atom(contact[0]).residue.resSeq # residue number of the first atom in the contact pair
        r2 = pdb.top.atom(contact[1]).residue.resSeq
        n1 = pdb.top.atom(contact[0]).residue.name
        n2 = pdb.top.atom(contact[1]).residue.name
        string = "{:d} {:s} {:} {:d} {:s} {:} {:3.5f}\n".format(r1, n1, "CA", r2, n2, "CA", filtered_dist1[i])
        fout.write(string)
    fout.write("\n") 
    # fout.write("\n") # empty line to separate contacts of pep1 from contacts of pep2
    
    for i, contact in enumerate(contact_pairs2):
        r1 = pdb.top.atom(contact[0]).residue.resSeq
        r2 = pdb.top.atom(contact[1]).residue.resSeq
        n1 = pdb.top.atom(contact[0]).residue.name
        n2 = pdb.top.atom(contact[1]).residue.name
        string = "{:d} {:s} {:} {:d} {:s} {:} {:3.5f}\n".format(r1, n1, "CA", r2, n2, "CA", filtered_dist2[i])
        fout.write(string)
        # fout.write("\n") 



with open ("pseudobonds_complex.pb", "w") as fout:
    for i, contact in enumerate(contact_pairs1):
        string = "/A:{:d}@CA /B:{:d}@CA\n".format(pdb.top.atom(contact[0]).residue.resSeq,
                                                  pdb.top.atom(contact[1]).residue.resSeq)
        fout.write(string)
    fout.write("\n")
    
    for i, contact in enumerate(contact_pairs2):
        string = "/A:{:d}@CA /B:{:d}@CA\n".format(pdb.top.atom(contact[0]).residue.resSeq,
                                                  pdb.top.atom(contact[1]).residue.resSeq)
        fout.write(string)