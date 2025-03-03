#!/bin/bash -l
#
#SBATCH --job-name=m1_new
#SBATCH --output=%x_%A_%a.o
#SBATCH --error=%x_%A_%a.e
#SBATCH --account=ppi
#SBATCH --partition=tier3
#SBATCH --time=24:00:00
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1g
#SBATCH --gpus=a100:15
#SBATCH --gpus-per-task=a100:1   # Cluster thinks 1 GPU per task, but MELD knows there are 2 tasks per GPU
#SBATCH --gpu-bind=single:1


spack env activate brini-x86_64-24082801

current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo "Run on $current_date at $current_time"

echo "Current Folder:"
pwd 
echo $PWD

echo "Spack Loaded Modules:"
spack find --loaded

# run the setup.py or restart the simulation, depending on whether simulation files exist. For a fresh run, remove Data/ and Logs/ and Logs/ folders.
# If a log file exists, the simulation continues. Otherwise, the simulation is set up by setup_MELD.py
if [ -e ./Logs/remd_000.log ]; then
  echo "Restarting MELD simulation"
  prepare_restart --prepare-run
else
  echo "Running setup_MELD.py..."
  python setup_MELD.py
fi

echo "Running the simulation"
srun launch_remd --debug # MELD-specifiic command to run the REMD simulation after looking for the info in Data/

