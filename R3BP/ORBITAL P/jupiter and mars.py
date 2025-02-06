import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from skyfield.api import load, Topos
import numpy as np

# Load planetary data
eph = load('de421.bsp')
ts = load.timescale()

# Define the planets
jupiter = eph['jupiter barycenter']
mars = eph['mars barycenter']
earth = eph['earth barycenter']

# Calculate positions over time
t = ts.utc(2023, 1, range(1, 366))  # Year 2023
jupiter_position = earth.at(t).observe(jupiter).position.km
mars_position = earth.at(t).observe(mars).position.km

# Extract x, y, and z coordinates
jupiter_x, jupiter_y, jupiter_z = jupiter_position
mars_x, mars_y, mars_z = mars_position

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot orbits
ax.plot(jupiter_x, jupiter_y, jupiter_z, label='Jupiter', color='orange')
ax.plot(mars_x, mars_y, mars_z, label='Mars', color='red')

# Plot the Sun
ax.scatter([0], [0], [0], color='yellow', marker='o', label='Sun')

# Customize the plot
ax.set_title('Orbits of Jupiter and Mars around the Sun (3D)')
ax.set_xlabel('X Position (km)')
ax.set_ylabel('Y Position (km)')
ax.set_zlabel('Z Position (km)')
ax.legend()

# Equal aspect ratio
ax.set_box_aspect([1, 1, 1])

# Show the plot
plt.show()
