import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Replace the following with your 20 sets of x, y, z coordinates
coordinates = [
    (2.271773623572096, 2.940633719465365, 0.2963756222591094),
    (2.57114740669013 , 2.5395544605502303, 0.12250936338189362),
    # Add more sets of coordinates here
]


# Extract x, y, and z coordinates from the list of tuples
x = [coord[0] for coord in coordinates]
y = [coord[1] for coord in coordinates]
z = [coord[2] for coord in coordinates]

# Create a 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the points
ax.scatter(x, y, z, c='b', marker='o')


# Set labels for the axes
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

# Show the plot
plt.show()
