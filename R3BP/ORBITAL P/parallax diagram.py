import matplotlib.pyplot as plt
import numpy as np

# Placeholder values for the Earth-Marion distance and Marion-Sun distance in AU.
# You should replace these with your actual calculated distances.
distance_to_marion_km = 75000000  # Replace with your actual value
distance_to_marion_au = 0.5       # Replace with your actual value
marion_sun_distance_au = 1.5      # Replace with your actual value
au_in_km = distance_to_marion_km / distance_to_marion_au

# Create figure and axis
fig, ax = plt.subplots()

# Plot the Earth
earth = plt.Circle((0, 0), 0.05, color='blue', label='Earth')
ax.add_patch(earth)

# Plot Marion's position
marion = plt.Circle((distance_to_marion_au, 0), 0.05, color='grey', label='Marion')
ax.add_patch(marion)

# Plot the Sun's position
sun = plt.Circle((marion_sun_distance_au, 0), 0.1, color='yellow', label='Sun')
ax.add_patch(sun)

# Draw a line representing the distance to Marion
plt.plot([0, distance_to_marion_au], [0, 0], color='black', linestyle='--')

# Annotate distances
plt.text(distance_to_marion_au/2, 0.1, f"{distance_to_marion_km} km", ha='center')
plt.text(distance_to_marion_au + marion_sun_distance_au/2, 0.1, f"{marion_sun_distance_au} AU", ha='center')

# Set the limits, labels, and title
ax.set_xlim(-1, 3)
ax.set_ylim(-2, 2)
ax.set_aspect('equal')
plt.xlabel('Astronomical Units (AU)')
plt.ylabel('Astronomical Units (AU)')
plt.title('Conceptual Plot for the Distance to Marion and Calculation of 1 AU')

# Add a legend
plt.legend()

# Show the plot
plt.show()
