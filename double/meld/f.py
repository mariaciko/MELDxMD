import numpy as np
import matplotlib.pyplot as plt

# Load the fup data from fup.dat
fup = np.loadtxt('fup.dat')

# Compute fdown
fdown = 1 - fup

# Generate the x-axis values (indices of the data points)
x_vals = np.arange(len(fup))

# Create a line plot
plt.figure(figsize=(12, 8))

# Plot fup and fdown
plt.plot(x_vals, fup, marker='o', linestyle='-', color='black', label='fup')
plt.plot(x_vals, fdown, marker='x', linestyle='--', color='r', label='fdown')

# Add diagonal reference lines
max_value = max(max(fup), max(fdown))
plt.plot([0, len(fup)-1], [0, max_value], linestyle='--', color='r', label='Diagonal 1')
plt.plot([0, len(fup)-1], [max_value, 0], linestyle='--', color='black', label='Diagonal 2')

# Labeling the plot
plt.xlabel('Replica Index')
plt.ylabel('Probability')
plt.title(r'Line Plot of $f_{up}$ and $f_{down}$: Probability of each runner moving up or down')
plt.legend()

# Save the plot to a file
plt.savefig('fup_fdown_line_plot.png')

# Display the plot
plt.show()
