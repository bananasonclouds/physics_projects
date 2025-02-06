import matplotlib.pyplot as plt
import numpy as np

#this data is for helio

def parse_ra_dec_individual(ra_str, dec_str):
    """Convert RA and Dec from 'hh mm ss.ss' and 'dd mm ss.ss' to decimal degrees."""
    ra_h, ra_m, ra_s = map(float, ra_str.split())
    dec_d, dec_m, dec_s = map(float, dec_str.split())

    ra_decimal = (ra_h + ra_m / 60 + ra_s / 3600) * 15
    dec_decimal = dec_d + dec_m / 60 + dec_s / 3600

    if dec_str.startswith('-'):
        dec_decimal = -dec_decimal

    return ra_decimal, dec_decimal

def parse_ra_dec(ra_dec_str):
    """Convert combined RA and Dec string to decimal degrees."""
    parts = ra_dec_str.split()
    ra_str = ' '.join(parts[:3])
    dec_str = ' '.join(parts[3:6])
    return parse_ra_dec_individual(ra_str, dec_str)

# Your observed data
observed_data_str = [
    "04 01 38.87 41 58 11.5", "04 00 59.79 41 57 26.8", "04 00 30.99 41 55 28.3",
    "03 59 33.62 41 51 00.8", "03 58 46.09 41 46 59.4", "03 55 40.98 41 29 29.1",
    "03 48 28.34 40 40 23.5", "03 46 46.57 40 27 22.5", "03 40 22.73 39 33 10.3",
    "03 33 36.02 38 25 55.7"
]

# Your reference data
reference_data_str = [
    "04 01 38.90 41 58 12.04", "04 00 59.690 41 57 26.64", "04 00 30.61 41 55 26.93",
    "03 59 33.32 41 50 58.89", "03 58 45.69 41 46 57.63", "03 55 40.78 41 29 27.60",
    "03 48 24.37 40 39 56.23", "03 46 40.16 40 26 34.12", "03 40 20.15 39 32 48.40",
    "03 33 23.14 38 23 42.16"
]

# Convert data to the required format
observed_data = [(str(i+1), data) for i, data in enumerate(observed_data_str)]
reference_data = [(str(i+1), data) for i, data in enumerate(reference_data_str)]

# Convert RA and Dec to decimal degrees
observed_data_decimal = [(time, *parse_ra_dec(ra_dec)) for time, ra_dec in observed_data]
reference_data_decimal = [(time, *parse_ra_dec(ra_dec)) for time, ra_dec in reference_data]

# Extracting time, RA, and Dec from the data
times_observed, ra_observed, dec_observed = zip(*observed_data_decimal)
times_reference, ra_reference, dec_reference = zip(*reference_data_decimal)

# Convert to numpy arrays for easier manipulation
times_observed = np.array(times_observed, dtype=float)
ra_observed = np.array(ra_observed)
dec_observed = np.array(dec_observed)

ra_reference = np.array(ra_reference)
dec_reference = np.array(dec_reference)

# Calculate residuals
ra_residuals = ra_observed - ra_reference
dec_residuals = dec_observed - dec_reference

# Plotting RA residuals against time
plt.figure(figsize=(9, 6))
plt.plot(times_observed, ra_residuals, marker='o', linestyle='-', color='blue')
plt.xlabel('Time')
plt.ylabel('RA Residuals (Degrees)')
plt.title('RA Residuals Over Time')
plt.grid(True)

# Plotting Dec residuals against time
plt.figure(figsize=(9,6))
plt.plot(times_observed, dec_residuals, marker='o', linestyle='-', color='red')
plt.xlabel('Time')
plt.ylabel('Dec Residuals (Degrees)')
plt.title('Dec Residuals Over Time')
plt.grid(True)

#plt.tight_layout()
plt.show()
