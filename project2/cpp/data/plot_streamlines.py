import numpy as np
import matplotlib.pyplot as plt

# Replace 'data.txt' with the path to your file
data = np.loadtxt('1000-128-stream')

# Extract columns
x = data[:, 0]  # x-coordinates
y = data[:, 1]  # y-coordinates
vorticity = data[:, 2]  # vorticity values
psi = data[:, 3]  # streamfunction values

# Determine the grid size (e.g., NxN grid)
N = int(np.sqrt(len(x)))

# Reshape the data into 2D arrays
x_grid = x.reshape(N, N)
y_grid = y.reshape(N, N)
psi_grid = psi.reshape(N, N)

# Create the plot
plt.figure(figsize=(8, 6))

# Contour plot for streamlines
plt.contour(x_grid, y_grid, psi_grid, levels=15, cmap='jet')  # Adjust 'levels' for desired detail
plt.colorbar(label='Stream Function (ψ)')

# Add labels and title
plt.title('Re=1000')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')  # Equal scaling for axes

# Show the plot
plt.show()
