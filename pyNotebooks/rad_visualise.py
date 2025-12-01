
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

# Load the data
data = np.load('poynting_data_ord_35.npz')
theta = data['theta']
phi = data['phi']
poynting_magnitudes = data['poynting_radial']

# Create a finer grid for interpolation
phi_grid, theta_grid = np.mgrid[0:2*np.pi:1000j, 0:np.pi:500j]  # Increased resolution

# Interpolate Poynting magnitudes onto the grid
points = np.vstack((phi, theta)).T
grid_values = griddata(points, poynting_magnitudes, (phi_grid, theta_grid), method='linear', fill_value=0.0)

# Convert grid to Cartesian coordinates for the sphere
x_grid = -np.sin(theta_grid) * np.cos(phi_grid)
y_grid = np.sin(theta_grid) * np.sin(phi_grid)
z_grid = np.cos(theta_grid)

# Debug: Check the range of interpolated values
print("Min grid_values:", np.nanmin(grid_values))
print("Max grid_values:", np.nanmax(grid_values))

# Create 3D surface plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(x_grid, y_grid, z_grid, facecolors=plt.cm.plasma(grid_values/np.nanmax(grid_values)), shade=False)
ax.set_title('Poynting Vector Magnitude on Unit Sphere')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Add a colorbar
mappable = plt.cm.ScalarMappable(cmap=plt.cm.plasma)
mappable.set_array(grid_values)
cbar = fig.colorbar(mappable, ax=ax, label='Poynting Vector Magnitude')

# Display the plot
plt.show()
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the data
data = np.load('poynting_data_ord_35.npz')
theta = data['theta']
phi = data['phi']
poynting_magnitudes = data['poynting_radial']

# Convert spherical to Cartesian coordinates
x = np.sin(theta) * np.cos(phi)
y = np.sin(theta) * np.sin(phi)
z = np.cos(theta)

# Create a 3D scatter plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(x, y, z, c=poynting_magnitudes, cmap='viridis', s=50)

# Customize the plot
ax.set_title('Color-coded Scatter Plot of Poynting Vector Data')
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

# Add a colorbar
cbar = fig.colorbar(scatter, ax=ax, label='Poynting Vector Magnitude')

# Display the plot
plt.show()
"""