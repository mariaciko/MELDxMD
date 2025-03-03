#!/usr/bin/env python
# encoding: utf-8
import numpy as np
import meld as meld
from meld import unit as u
from meld.system.scalers import LinearRamp
from meld.remd import ladder, adaptor, leader
from meld import comm, vault, system, parse
import glob as glob
from openmm import *
from openmm import app as app
from openmm import unit as u
from openmm import LangevinIntegrator, MonteCarloBarostat
import mdtraj as md

N_REPLICAS = 30
N_STEPS = 150000
BLOCK_SIZE = 100

# Load the system
templates = glob.glob('TEMPLATES/*.pdb')
print(templates)
traj = md.load('TEMPLATES/pep_shifted.pdb')
print(f'Number of frames in trajectory: {traj.n_frames}')

# Set up the system
build_options = meld.AmberOptions(
    forcefield = 'ff14sbside',
    implicit_solvent_model='gbNeck2',
    use_bigger_timestep = True, 
    cutoff= 1.8*u.nanometer)
p = meld.AmberSubSystemFromPdbFile(templates[0])
b = meld.AmberSystemBuilder(build_options)
s = b.build_system([p]).finalize()
s.temperature_scaler = meld.GeometricTemperatureScaler(0, 0.4, 300. *u.kelvin, 500.* u.kelvin)
print("system ok")

# specify chains
protein_chain = traj.topology.select('chainid == 0')
peptide1_chain = traj.topology.select('chainid == 1')
peptide2_chain = traj.topology.select('chainid == 2')
CA_peptide1 = traj.topology.select('chainid == 1 and name == CA')
CA_peptide2 = traj.topology.select('chainid == 2 and name == CA')
print(f'Length of peptide1 chain: {len(CA_peptide1)}\nLength of peptide2 chain: {len(CA_peptide2)}')
CA1_central_index = CA_peptide1[6] 
CA2_central_index = CA_peptide2[len(CA_peptide2) // 2]
print(f'CA1 central index: {CA1_central_index}\nCA2 central index: {CA2_central_index}')

n_res_protein = len(set(traj.topology.atom(i).residue.index for i in protein_chain))
n_res_peptide1 = len(set(traj.topology.atom(i).residue.index for i in peptide1_chain))
n_res_peptide2 = len(set(traj.topology.atom(i).residue.index for i in peptide2_chain))
n_res = n_res_protein + n_res_peptide1 + n_res_peptide2
print(f'Number of protein residues: {n_res_protein}\nNumber of peptide1 residues: {n_res_peptide1}\nNumber of peptide2 residues: {n_res_peptide2}\nTotal number of residues: {n_res}')

rg_protein = md.compute_rg(traj.atom_slice(protein_chain)) * 0.1  # convert Å to nm
rg_peptide1 = md.compute_rg(traj.atom_slice(peptide1_chain)) * 0.1
rg_peptide2 = md.compute_rg(traj.atom_slice(peptide2_chain)) * 0.1
print(f'Radius of gyration for protein: {rg_protein} nm\nRadius of gyration for peptide1: {rg_peptide1} nm \nRadius of gyration for peptide2: {rg_peptide2} nm')

# sphere radius
cutoff_distance = 1.8 # nm
sphere_radius = 7.0*u.nanometer
protein_com = md.compute_center_of_mass(traj.atom_slice(protein_chain)) # nm
sphere_center = protein_com  # align sphere center with protein COM
print(f"Radius of the sphere: {sphere_radius}\nProtein center of mass: {protein_com}\nUpdated sphere center: {sphere_center}")
peptide_coms = [
    md.compute_center_of_mass(traj.atom_slice(peptide1_chain)), # nm
    md.compute_center_of_mass(traj.atom_slice(peptide2_chain))] # nm
print(f"Peptide1 com: {peptide_coms[0]} \nPeptide2 com: {peptide_coms[1]}")

topology_filename = 'topology_setup.pdb'
traj.save(topology_filename)
print(f'Topology file saved as {topology_filename}')


# Distance restraint function for peptide
def get_dist_restraints(filename, s, scaler): # file has atom1, atom2, distance
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    # Do an if test: if the information is there, read it and create a restraint. If not, make a group.
    for line in lines:
        if not line: # if line is empty, make a group (enforce 1 of 2 groups). Otherwise, keep adding to the group
            if rest_group:
                dists.append(s.restraints.create_restraint_group(rest_group, 1))
                rest_group = []
        else:
            cols = line.split() 
            i = int(cols[0])
            name_i = "CA"
            j = int(cols[3])
            name_j = "CA"
            dist = float(cols[6])
            atom_1_index = s.index.atom(i, name_i)
            atom_2_index = s.index.atom(j, name_j)
            rest = s.restraints.create_restraint('distance', scaler, LinearRamp(0.0,100.0,0.0,1.0),   # scaler turns off the springs at α = 0.4, LinearRamp turns on the springs at the beginning of sim
              r1=0.0*u.nanometer, r2=0.0*u.nanometer, r3=dist*u.nanometer, # r2 and r3 is the width where potential is 0
              r4=(dist+0.2)*u.nanometer, k=250*u.kilojoules_per_mole/(u.nanometer*u.nanometer),
              atom1=atom_1_index, atom2=atom_2_index)
            rest_group.append(rest)
    return dists

# Cartesian restraint function to keep protein rigid
def make_cartesian_collections(s, scaler, residues, delta=0.2*u.nanometers, k=250.):
    cart = []
    backbone = ['CA'] 
    for i in residues:
        for b in backbone:
            atom_index = s.index.atom(i,b)
            x,y,z = (s.template_coordinates[atom_index]/10)*u.nanometers
            rest = s.restraints.create_restraint('cartesian', scaler, LinearRamp(0.0, 15.0, 0.0, 1.0), atom_index=atom_index, 
                x=x, y=y, z=z, delta=delta,force_const=k*u.kilojoules_per_mole/(u.nanometer*u.nanometer))
            # xyz must be fixed; delta = width of flat bottom potential (how much we allow the Cartesian restraint to fluctuate)
            cart.append(rest)*u.nanometer
    return cart

# Confinement restraint function to keep peptide within sphere
def make_sphere_distance_restraints(s, scaler, peptide_residues, sphere_radius, force_const):
    dist = []
    backbone = ['CA']
    for i in peptide_residues:
        for b in backbone:
            atom_index = s.index.atom(i, b) 
            rest = meld.system.restraints.ConfinementRestraint(s, scaler, LinearRamp(0.0,100.0,0.0,1.0),
                atom_index=atom_index, radius=sphere_radius,
                force_const=force_const*u.kilojoules_per_mole/(u.nanometer*u.nanometer))
            dist.append(rest)
    return dist

# Exclusion restraint function to keep peptides from overlapping
def pep_pep_exclusion(s, scaler, CA1_central_index, CA2_central_index, force_const):
    rest = s.restraints.create_restraint('distance', scaler, LinearRamp(0.0,100.0,0.0,1.0),
        r1=0.0*u.nanometer, r2=3.0*u.nanometer, r3=100*u.nanometer, r4=101*u.nanometer, 
        k=force_const*u.kilojoules_per_mole/(u.nanometer*u.nanometer),
        atom1=CA1_central_index, atom2=CA2_central_index)
    return rest

# Keep one peptide in reference state at low alpha at any given time
def prot_pep_ref(s, scaler, peptide1_residues, peptide2_residues, protein_residues, force_const):
    rest_list = []
    for peptide_com in peptide_coms:
        rest = s.restraints.create_restraint('distance', scaler, LinearRamp(0.0, 100.0, 0.0, 1.0),
                                          r1=4.0 * u.nanometer, r2=5.0 * u.nanometer, r3=7.0 * u.nanometer, 
                                          r4=8.0 * u.nanometer, k=force_const * u.kilojoules_per_mole / (u.nanometer * u.nanometer),
                                          atom1=protein_com, atom2=peptide_com)
        rest_list.append(rest)

    # Create a MELD collection to enforce the restraint with the lowest energy
    ref_collection = s.restraints.create_collection(rest_list, 'conditional_lowest_energy', 
                                                    condition=lambda r1, r2: r1 if r1.energy < r2.energy else r2)
    return ref_collection
    


def setup_system():
    # cartesian restraint to keep protein fixed
    prot_scaler = s.restraints.create_scaler('constant')
    prot_cart_rest = make_cartesian_collections(s, prot_scaler, np.arange(1, n_res_protein+1))
    s.restraints.add_as_always_active_list(prot_cart_rest)
    print('Protein restraints (prot_cart_rest):', len(prot_cart_rest))
    
    # distance restraints to keep peptide bound to protein
    prot_pep_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=0.9, factor=2.0)
    prot_pep_rest = get_dist_restraints('protein_pep_all.dat',s,scaler=prot_pep_scaler)
    s.restraints.add_selectively_active_collection(prot_pep_rest, 1)
    print('Protein-peptide restraints (prot_pep_rest):', len(prot_pep_rest))

    # sphere restraint
    sphere_scaler = s.restraints.create_scaler('constant')
    pep_res = list(range(n_res_protein, n_res))
    sphere_rest = make_sphere_distance_restraints(s, sphere_scaler, pep_res, sphere_radius, force_const=250)
    s.restraints.add_as_always_active_list(sphere_rest)
    print('Sphere peptide restraints (sphere_rest):', len(sphere_rest))

    # peptide-peptide exclusion restraint
    exclusion_scaler = s.restraints.create_scaler('constant')
    # pep1_res = list(range(n_res_protein, n_res_protein+n_res_peptide1))
    # pep2_res = list(range(n_res_protein+n_res_peptide1, n_res))
    pep_excl_rest = pep_pep_exclusion(s, exclusion_scaler, CA1_central_index, CA2_central_index, force_const=250)
    s.restraints.add_as_always_active_list(pep_excl_rest)
    print('Peptide-peptide restraints (pep_excl_rest):', len(pep_excl_rest))

    # reference state restraint
    ref_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.0, alpha_max=0.6, factor=2.0)
    ref_rest = prot_pep_ref(s, ref_scaler, peptide1_chain, peptide2_chain, protein_chain, force_const=250)
    s.restraints.add_as_always_active_list(ref_rest)
    print('Reference state restraints (ref_rest):', len(ref_rest))

    options = meld.RunOptions(
        timesteps = 11111,
        minimize_steps = 20000)

    remd = meld.setup_replica_exchange(s, n_replicas=N_REPLICAS, n_steps=N_STEPS)
    meld.setup_data_store(s, options, remd)
    print('store ok')

setup_system()