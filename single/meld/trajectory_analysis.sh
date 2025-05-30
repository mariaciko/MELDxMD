#!/bin/bash -l
#
#SBATCH --job-name=analysis
#SBATCH --output=%x_%A_%a.o
#SBATCH --error=%x_%A_%a.e
#SBATCH --account=ppi
#SBATCH --partition=tier3
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1g

spack env activate brini-x86_64-24082801

echo "Current Folder:"
pwd 
echo $PWD

echo "Spack Loaded Modules:"
spack find --loaded

# Run the required Python scripts for analysis
echo "Extracting trajectory for replica 0"
extract_trajectory extract_traj_dcd --replica 0 traj.dcd
echo "Extracting trajectory for replica 29"
extract_trajectory extract_traj_dcd --replica 29 traj29.dcd

echo "Extracting trajectory end for replica 0"
extract_trajectory extract_traj --end 2 --replica 0 topology.pdb
echo "Extracting trajectory end for replica 29"
extract_trajectory extract_traj --end 2 --replica 29 topology29.pdb

echo "Extracting trace data"
analyze_remd extract_trace trace.dat
python analyze_remd.py visualize_trace
python analyze_remd.py visualize_accept
echo "Visualized trace and accept"

echo "Extracting fup data"
analyze_remd extract_fup fup.dat
python analyze_remd.py visualize_fup
echo "Visualized fup"

echo "Extracting alpha data"
analyze_remd extract_alpha alpha.dat
python analyze_remd.py visualize_alpha
echo "Visualized alpha"


echo "Trajectory analysis complete"


echo "Generating plots..."
python distance.py
python RMSD.py
python roundtrips.py
echo "Plots generated"
