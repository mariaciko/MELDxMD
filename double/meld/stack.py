#!/usr/bin/env python
# coding: utf-8
import numpy as np
import mdtraj as md

pdb1 = md.load('../native.pdb')
pdb2 = md.load('../peptide2.pdb')

protein = pdb1.atom_slice(pdb1.top.select("chainid == 0"))
peptide1 = pdb1.atom_slice(pdb1.top.select("chainid == 1"))
peptide2 = pdb2.atom_slice(pdb2.top.select("chainid == 0"))

# Stack the two trajectories
stacked_pdb = pdb1.stack(peptide2)

stacked_pdb.save_pdb('stacked.pdb')

print("Stacked trajectory saved as 'stacked.pdb'")
