import matplotlib.pyplot as plt
import numpy as np

def plot_orbit(semi_major_axis, color, label):
    theta = np.linspace(0, 2*np.pi, 100)
    x = semi_major_axis * np.cos(theta)
    y = semi_major_axis * np.sin(theta)
    plt.plot(x, y, color=color, label=label)

def plot_main_belt(inner, outer, color):
    theta = np.linspace(0, 2*np.pi, 100)
    for r in np.linspace(inner, outer, 30):  # Creates the annular effect
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        plt.plot(x, y, color=color, alpha=0.3)

# Define semi-major axes (in arbitrary units)
earth_orbit = 1.0  # Assuming Earth's orbit as 1 AU
mars_orbit = 1.5
jupiter_orbit = 5.2
asteroid_belt_inner = 2.1
asteroid_belt_outer = 3.3

# Plot the orbits
plt.figure(figsize=(10, 10))
plot_orbit(earth_orbit, 'blue', 'Earth Orbit')
plot_orbit(mars_orbit, 'red', 'Mars Orbit')
plot_orbit(jupiter_orbit, 'orange', 'Jupiter Orbit')
plot_main_belt(asteroid_belt_inner, asteroid_belt_outer, 'grey')

# Plot the planets
plt.scatter([earth_orbit], [0], color='blue', s=40)  # Earth
plt.scatter([mars_orbit], [0], color='red', s=50)  # Mars
plt.scatter([jupiter_orbit], [0], color='orange', s=80)  # Jupiter

# Plot the Sun
plt.scatter([0], [0], color='yellow', s=100, label='Sun')

# Manually add a legend entry for the Main Belt
plt.plot([], [], color='grey', alpha=0.3, label='Main Belt')

# Set up plot labels and legend
#plt.title("Solar System: Sun, Earth, Mars, Jupiter, and Main Belt")
plt.xlabel("Distance from Sun (AU)")
plt.ylabel("Distance from Sun (AU)")
plt.legend(loc='upper right')
plt.grid(True)
plt.axis('equal')

plt.show()
