## MELDxMD protein-peptide binding

We want to predict binding modes using MELD-accelerated Molecular Dynamics, where MELD stands for Modeling Employing Limited Data. This approach accelerates conventional MD methods by integrating external information in a Bayesian framework. A common issue with physics-based approaches like MD are the kinetic traps conformations fall into in the potential energy surface, resulting in a process that is slow and computationally expensive. 

This project focuses on predicting the binding poses for P53-derived ligands to the MDM-2, X protein using MELDxMD via two methods: 

(i) predict binding in a competitive environment;  
(ii) predict binding by solving the problem of steric hindrance

For (i), we follow the study of _Morrone, J. A., Perez, A., MacCallum, J. & Dill, K. A. Computed Binding of Peptides to Proteins with MELD-Accelerated Molecular Dynamics. J. Chem. Theory Comput. 13, 870–876 (2017)_. Method (ii) was applied to the nine complexes outilned in Morrone's study. This repository contains the preparatory files for the native structure 1YCR (MDM2-P53) and the reference peptide (P6W) for method (i) and (ii), in [single](https://github.com/mariaciko/Thesis/tree/main/single) and [double](https://github.com/mariaciko/Thesis/tree/main/double), respectively. 

Prior to running the simulations on MELD, we test for convergence running a single-trajectory 250-million step MD simulation for each complex in OpenMM and verify energy stabilization by consulting RMSD distributions of the receptor, the ligand, and the ligand whilst holding the receptor aligned, checking if (upon binding) the protein unfolded, the ligand changed shape, and if unbinding occurred. All of the systems passed this check and so all were cleared to be passed on to MELD.

### MELD protocol
MELD relies on Bayesian statistics to narrow its conformational search to only desired states. It implements "smart springs" or restraints that are not all active at the same time though various scalers (e.g. protein cartesian restraints, protein-peptide binding, containment restraints, etc).  The idea behind "smart" springs is that they rely on a flat bottom potential and Bayesian statistics.
1. Flat bottom potential
Take a quadratic spring. Cut it in half. Distance the halves and bring them down to the horizontal axis. Now, you have a flat bottom region which represents the zero-potential. This means that points falling in this flat region will behave as though they are at equilibrium. Equilibrium matters because that's where you compute thermodynamic properties from. If we were to work with the original quadratic spring, equilibrium only occurs at one point: the minimum of the curve. 

2. Bayesian statistics
Bayes' Theorem states that the springs which are active are the least violated ones (i.e. low-energy conformations).
<img width="380" alt="image" src="https://github.com/user-attachments/assets/8790cd18-7729-4617-a729-69cc38318a11" />
In MELD, the probability of sampling a configuration given a set of springs is defined as the probability of sampling a configuration in MD (prior) multiplied by the probability of having a spring active (likelihood).
<img width="335" alt="Screenshot 2025-05-10 at 12 18 25" src="https://github.com/user-attachments/assets/875ba9d0-0ba8-43f5-b683-acb1cdf5dd11" />


MELD also uses Hamiltonian and Temperature Replica Exchange Molecular Dynamics (H, T-REMD) to maximize sampling efficiently. While the landscape anneals down to sample narrower distibutions closer to the native state, replica exchange uses copies of the system, or replicates, to explore conformations that MELD would otherwise penalize.

