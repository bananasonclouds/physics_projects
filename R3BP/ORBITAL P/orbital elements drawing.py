import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import table

# Data
data = {
    "Asteroid": ["Helio", "Aneas", "Kallisto", "Deborah", "Marion", "Eurykleia", "Themis"],
    "Orbital Period (Years)": [5.75, 11.96, 4.38, 4.74, 5.31, 4.90, 5.59],
    "Diameter (km)": [109.60, 118.0, 48.60, 60.10, 105.90, 93.10, 198.0],
    "Rotation Period (Hours)": [234840.0, 8.71, 19.49, 29.37, 13.55, 16.52, 8.37],
    "Spectral Type (Tholen/SMASS)": ["B", "D", "S", "B", "XC", "C", "C"],
    "Hazard Level": ["Not Hazardous", "Not Hazardous", "Not Hazardous", "Not Hazardous", "Not Hazardous", "Not Hazardous", "Not Hazardous"]
}

# Create DataFrame
df = pd.DataFrame(data)

# Plotting
fig, ax = plt.subplots(figsize=(12, 4))  # Adjust size as needed
ax.axis('tight')
ax.axis('off')
the_table = table(ax, df, loc='center')

#plt.title('Asteroid Information Table')

plt.show()
