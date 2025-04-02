# Thesis
## MELDxMD protein-peptide binding

We want to predict binding modes using MELD-accelerated Molecular Dynamics, where MELD stands for Modeling Employing Limited Data. This approach accelerates conventional MD methods by integrating external information in a Bayesian framework. A common issue with physics-based approaches like MD are the kinetic traps conformations fall into in the potential energy surface, resulting in a process that is slow and computationally expensive. 

This project focuses on predicting the binding poses for P53-derived ligands to the MDM-2, X protein using MELDxMD via two methods: 

(i) predict binding in a competitive environment;  
(ii) predict binding by solving the problem of steric hindrance

For (i), we follow the study of __Morrone, J. A., Perez, A., MacCallum, J. & Dill, K. A. Computed Binding of Peptides to Proteins with MELD-Accelerated Molecular Dynamics. J. Chem. Theory Comput. 13, 870–876 (2017)__. This repository contains the preparatory files for the native structure 1YCR (MDM2-P53) and the reference peptide (P6W) for method (i) and (ii), in [single](https://github.com/mariaciko/Thesis/tree/main/single) and [double](https://github.com/mariaciko/Thesis/tree/main/double), respectively.  

Prior to running the simulations on MELD, we test for convergence running a 250-million step MD simulation for each complex in OpenMM and verify energy stabilization by consulting RMSD distributions of the receptor, the ligand, and the ligand whilst holding the receptor aligned. All of the systems passed this check and so all were cleared to be passed on to MELD.

### MELD protocol
MELD relies on Bayesian statistics to narrow its conformational search to only desired states (i.e. low-energy conformations). It implements "smart springs" or restraints that are not active indefinitely though various scalers (e.g. protein cartesian restraints, protein-peptide binding, containment restraints, etc). MELD also used Hamiltonian and Temperature Replica Exchange Molecular Dynamics (H, T-REMD) to maximize sampling efficiently. While the landscape anneals down to sample narrower distibutions closer to the native state, replica exchange uses copies of the system, or replicates, to explore conformations that MELD would otherwise penalize.

