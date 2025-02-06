import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

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
    "4 01 38.90 41 58 12.04", "04 00 59.690 41 57 26.64", "04 00 30.61 41 55 26.93",
    "03 59 33.32 41 50 58.89", "03 58 45.69 41 46 57.63", "03 55 40.78 41 29 27.60",
    "03 48 24.37 40 39 56.23", "03 46 40.16 40 26 34.12", "03 40 20.15 39 32 48.40",
    "03 33 23.14 38 23 42.16"
]

# Convert data to the required format
observed_data = [parse_ra_dec(data) for data in observed_data_str]
reference_data = [parse_ra_dec(data) for data in reference_data_str]

# Extracting RA and Dec from the data
ra_observed, dec_observed = zip(*observed_data)
ra_reference, dec_reference = zip(*reference_data)

# Calculate changes in RA and Dec
delta_ra = [o - r for o, r in zip(ra_observed, ra_reference)]
delta_dec = [o - r for o, r in zip(dec_observed, dec_reference)]

# Adjust RA change by the cosine of average Declination
avg_dec = [(o + r) / 2 for o, r in zip(dec_observed, dec_reference)]
delta_ra_cos_dec = [d_ra * np.cos(np.radians(d_dec)) for d_ra, d_dec in zip(delta_ra, avg_dec)]

# Assuming uncertainties in RA and Dec
# Replace these with your actual uncertainties
sigma_ra = 0.1  # RA uncertainty in arcsec
sigma_dec = 0.1  # Dec uncertainty in arcsec

# Convert uncertainties from arcsec to degrees
sigma_ra /= 3600
sigma_dec /= 3600

# Adjust RA uncertainty by the average cosine of the Declination
sigma_ra_cos_dec = sigma_ra * np.mean([np.cos(np.radians(d)) for d in dec_observed])

# Plotting the Astrometric Position Comparison Plot
plt.figure(figsize=(8, 8))
plt.scatter(delta_ra_cos_dec, delta_dec, color='black', marker='.')

# Adding error ellipses for 1, 2, and 3 standard deviations
for nsigma in [1, 2, 3]:
    ellipse = patches.Ellipse((0, 0), 2 * nsigma * sigma_ra_cos_dec, 2 * nsigma * sigma_dec,
                              edgecolor='grey', facecolor='none', linestyle='dashed', linewidth=1, label=f'{nsigma}σ')
    plt.gca().add_patch(ellipse)

# Plotting a red dot at the mean position
plt.plot(np.mean(delta_ra_cos_dec), np.mean(delta_dec), 'ro')

plt.xlabel('Delta in RA*cos(Dec) [arcsec]')
plt.ylabel('Delta in Dec [arcsec]')
plt.title('Astrometric Position Comparison Plot')
plt.axhline(0, color='grey', linestyle='dashed')
plt.axvline(0, color='grey', linestyle='dashed')
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', 'box')
plt.show()
