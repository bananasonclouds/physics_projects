import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Given system parameters
G = 6.6726e-11  # m^3 kg^-1 s^-2
M_sun = 1.99e30  # kg
m1 = 6 * M_sun  # Star 1 mass
m2 = 2 * M_sun  # Star 2 mass
a_min = 7.75e11 # Minimum stable separation, 1 AU for simplicity, adjust as needed
a_min2= 8e11
a_max= 15e11
a_max1=20e11
# Define a grid for the plot
x = np.linspace(-25e11, 25e11, 400)
y = np.linspace(-25e11, 25e11, 400)
X, Y = np.meshgrid(x, y)

# Compute separation from the COM for each point
# Assuming COM is at (0, 0) for simplicity
R = np.sqrt(X**2 + Y**2)

# Normalize the separation distances to range [0, 1] for the colormap,
# where 1 indicates at or beyond a_min (stable), and 0 indicates far from a_min (unstable)
normalized_separation1 = np.clip((R - a_min) / (np.max(R) - a_min), 0, 1)


# Plot
#sns.color_palette("icefire", as_cmap=True)
plt.figure(figsize=(10, 8))
contour = plt.contourf(X, Y, normalized_separation1, levels=100,cmap=sns.color_palette("Spectral", as_cmap=True))
plt.plot(0, 0, 'yo', markersize=10, label='Center of Mass',color='black')  # Mark the COM
cbar = plt.colorbar(contour)
cbar.set_label('Normalized Separation Distance', fontsize=12)
plt.axhline(0, color='white', lw=1,linestyle='-')  # Add a line to represent the separation axis
plt.axvline(a_min, color='black', linestyle='--', label='$a_{min}$ Threshold for μ=0.25')
plt.axvline(a_min2, color='white', linestyle='--', label='$a_{min}$ Threshold for μ=0.5')
plt.axvline(a_max, color='black', linestyle='--', label='$a_{max}$ Threshold for μ=0.25')
plt.axvline(a_max1, color='white', linestyle='--', label='$a_{max}$ Threshold for μ=0.5')
plt.xlabel('X (m)',fontsize=16)
plt.ylabel('Y (m)',fontsize=16)
#plt.title('Stability Gradient in a Binary Star System')
plt.legend(loc='upper left', fontsize='large', markerscale=1.5)
#plt.grid(True)
plt.show()
