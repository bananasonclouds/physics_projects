import matplotlib.pyplot as plt
import numpy as np

# Example data for 'helio'
# Replace these example values with your actual data
observed_values = np.array([177.25, 145.8655, 55.74])  # Observed values
jpl_values = np.array([178.054, 146.18, 55.48])       # JPL values
errors = np.array([0.5, 0.3, 0.2])                    # Errors for observed values

# Number of observations
n = len(observed_values)
observation_numbers = np.arange(1, n + 1)

# Plot observed values with error bars
plt.errorbar(observation_numbers, observed_values, yerr=errors, fmt='o', label='Observed', capsize=5)

# Plot JPL values
plt.errorbar(observation_numbers, jpl_values, yerr=np.zeros_like(jpl_values), fmt='^', label='JPL', capsize=5)

# Add legend, labels, title, and grid
plt.legend()
plt.xlabel('Observation Number')
plt.ylabel('Argument of Periapsis')
plt.title('Comparison of Argument of Periapsis for helio')
plt.grid(True)

# Show plot
plt.show()
