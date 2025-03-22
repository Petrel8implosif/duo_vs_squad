import numpy as np
import matplotlib.pyplot as plt

# Load data from the text file
data = np.loadtxt('time_and_size.txt')  # Replace 'data.txt' with the actual file name
matrix_size = data[:, 0]
execution_time = data[:, 1]

# Create a log-log plot
plt.figure()
plt.loglog(matrix_size, execution_time, marker='o', linestyle='-', color='b', label='Execution Time')

# Generate O(n^3) and O(n^2) lines for comparison
n_squared = matrix_size**2
n_cubed = matrix_size**3

# Normalize the O(n^2) and O(n^3) lines for better visualization
n_squared_normalized = n_squared * (execution_time[0] / n_squared[0])
n_cubed_normalized = n_cubed * (execution_time[0] / n_cubed[0])

# Plot the O(n^2) and O(n^3) lines
plt.loglog(matrix_size, n_squared_normalized, linestyle='--', color='g', label='O(n^2)')
plt.loglog(matrix_size, n_cubed_normalized, linestyle='--', color='r', label='O(n^3)')

# Add labels, title, and legend
plt.xlabel('Matrix Size')
plt.ylabel('Execution Time (s)')
plt.title('Log-Log Plot of Matrix Size vs Execution Time')
plt.legend()

# Show the plot
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()