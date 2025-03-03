#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def follow_replica_ladder(trace_file, runner_id, output_plot):
    """
    Follows a simulation from the ground replica (runner_id=0) as it climbs the replica ladder.

    Parameters:
    trace_file (str): Path to the trace.dat file from the MELD simulation.
    runner_id (int): Runner ID to track during the simulation.

    Returns:
    None
    """
    # Load the trace.dat file
    trace = np.loadtxt(trace_file)
    
    # Ensure the trace file is transposed if needed
    if trace.shape[1] < trace.shape[0]:
        trace = trace.T
    
    # Ensure the runner_id is valid
    if runner_id < 0 or runner_id >= trace.shape[1]:
        print(f"Invalid runner_id: {runner_id}. Must be between 0 and {trace.shape[1] - 1}.")
        return
    
    # Use the raw steps as the x-axis
    steps = np.arange(0, len(trace[1]))
    print(trace.shape) 
    replica_positions = trace[runner_id, :]
    print(replica_positions.shape)


    # Create the step plot
    plt.figure(figsize=(10, 6))
    plt.step(steps, replica_positions, where='post', label=f'Runner {runner_id}', color='b')
    plt.xlabel('Steps')  # Change label to Steps
    plt.ylabel('Runner Position')
    plt.title(f'Runner {runner_id} Movement over Time')
    plt.xticks(np.arange(0, steps[-1] + 1, 15))  # Set ticks every 10 ps
    plt.grid(True)

    # Save the plot
    plt.savefig(output_plot)
    print(f"Plot saved to {output_plot}")


# Example usage if running as a standalone script:
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: ./follow_replica_ladder.py <trace_file> <runner_id> <output_plot>")
        sys.exit(1)

trace_file = sys.argv[1]
runner_id = int(sys.argv[2])
output_plot = sys.argv[3]

follow_replica_ladder(trace_file, runner_id, output_plot)
