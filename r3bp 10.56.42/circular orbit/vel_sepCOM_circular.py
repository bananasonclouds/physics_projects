import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.6726 * pow(10, -11)  # m^3 kg^-1 s^-2
M_sun = 1.99 * pow(10, 30)  # in kg

# Masses of the two stars and the planet (in kg)
m1 = 4.95 * M_sun  # Star 1
m2 = 2.20 * M_sun  # Star 2
m3 = 2 * pow(10, 27)  # Jupiter mass

miu = m2 / (m1 + m2)  # mass ratio

# Acceleration function including the planet
def acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3):
    # Calculate distances between bodies
    r12 = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    r13 = np.sqrt((x3 - x1) ** 2 + (y3 - y1) ** 2)
    r23 = np.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)

    # Calculate accelerations
    ax1 = G * m2 * (x2 - x1) / r12 ** 3 + G * m3 * (x3 - x1) / r13 ** 3
    ay1 = G * m2 * (y2 - y1) / r12 ** 3 + G * m3 * (y3 - y1) / r13 ** 3

    ax2 = G * m1 * (x1 - x2) / r12 ** 3 + G * m3 * (x3 - x2) / r23 ** 3
    ay2 = G * m1 * (y1 - y2) / r12 ** 3 + G * m3 * (y3 - y2) / r23 ** 3

    ax3 = G * m1 * (x1 - x3) / r13 ** 3 + G * m2 * (x2 - x3) / r23 ** 3
    ay3 = G * m1 * (y1 - y3) / r13 ** 3 + G * m2 * (y2 - y3) / r23 ** 3

    return ax1, ay1, ax2, ay2, ax3, ay3

# Runge-Kutta function including the planet
def runge_kutta(m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt):
    ax1, ay1, ax2, ay2, ax3, ay3 = acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3)

    # First set of Runge-Kutta
    kx1_1 = vx1 * dt
    ky1_1 = vy1 * dt
    kvx1_1 = ax1 * dt
    kvy1_1 = ay1 * dt
    kx2_1 = vx2 * dt
    ky2_1 = vy2 * dt
    kvx2_1 = ax2 * dt
    kvy2_1 = ay2 * dt
    kx3_1 = vx3 * dt
    ky3_1 = vy3 * dt
    kvx3_1 = ax3 * dt
    kvy3_1 = ay3 * dt

    # Intermediate positions for the second set
    x1_mid = x1 + kx1_1 / 2
    y1_mid = y1 + ky1_1 / 2
    x2_mid = x2 + kx2_1 / 2
    y2_mid = y2 + ky2_1 / 2
    x3_mid = x3 + kx3_1 / 2
    y3_mid = y3 + ky3_1 / 2

    # Second set of Runge-Kutta
    ax1_mid, ay1_mid, ax2_mid, ay2_mid, ax3_mid, ay3_mid = acceleration(
        m1, m2, m3, x1_mid, y1_mid, x2_mid, y2_mid, x3_mid, y3_mid
    )

    kx1_2 = (vx1 + kvx1_1 / 2) * dt
    ky1_2 = (vy1 + kvy1_1 / 2) * dt
    kvx1_2 = ax1_mid * dt
    kvy1_2 = ay1_mid * dt
    kx2_2 = (vx2 + kvx2_1 / 2) * dt
    ky2_2 = (vy2 + kvy2_1 / 2) * dt
    kvx2_2 = ax2_mid * dt
    kvy2_2 = ay2_mid * dt
    kx3_2 = (vx3 + kvx3_1 / 2) * dt
    ky3_2 = (vy3 + kvy3_1 / 2) * dt
    kvx3_2 = ax3_mid * dt
    kvy3_2 = ay3_mid * dt

    # Intermediate positions for the third set
    x1_mid = x1 + kx1_2 / 2
    y1_mid = y1 + ky1_2 / 2
    x2_mid = x2 + kx2_2 / 2
    y2_mid = y2 + ky2_2 / 2
    x3_mid = x3 + kx3_2 / 2
    y3_mid = y3 + ky3_2 / 2

    # Third set of Runge-Kutta
    ax1_mid, ay1_mid, ax2_mid, ay2_mid, ax3_mid, ay3_mid = acceleration(
        m1, m2, m3, x1_mid, y1_mid, x2_mid, y2_mid, x3_mid, y3_mid
    )

    kx1_3 = (vx1 + kvx1_2 / 2) * dt
    ky1_3 = (vy1 + kvy1_2 / 2) * dt
    kvx1_3 = ax1_mid * dt
    kvy1_3 = ay1_mid * dt
    kx2_3 = (vx2 + kvx2_2 / 2) * dt
    ky2_3 = (vy2 + kvy2_2 / 2) * dt
    kvx2_3 = ax2_mid * dt
    kvy2_3 = ay2_mid * dt
    kx3_3 = (vx3 + kvx3_2 / 2) * dt
    ky3_3 = (vy3 + kvy3_2 / 2) * dt
    kvx3_3 = ax3_mid * dt
    kvy3_3 = ay3_mid * dt

    # Intermediate positions for the fourth set
    x1_end = x1 + kx1_3
    y1_end = y1 + ky1_3
    x2_end = x2 + kx2_3
    y2_end = y2 + ky2_3
    x3_end = x3 + kx3_3
    y3_end = y3 + ky3_3

    # Fourth set of Runge-Kutta
    ax1_end, ay1_end, ax2_end, ay2_end, ax3_end, ay3_end = acceleration(
        m1, m2, m3, x1_end, y1_end, x2_end, y2_end, x3_end, y3_end
    )

    kx1_4 = (vx1 + kvx1_3) * dt
    ky1_4 = (vy1 + kvy1_3) * dt
    kvx1_4 = ax1_end * dt
    kvy1_4 = ay1_end * dt
    kx2_4 = (vx2 + kvx2_3) * dt
    ky2_4 = (vy2 + kvy2_3) * dt
    kvx2_4 = ax2_end * dt
    kvy2_4 = ay2_end * dt
    kx3_4 = (vx3 + kvx3_3) * dt
    ky3_4 = (vy3 + kvy3_3) * dt
    kvx3_4 = ax3_end * dt
    kvy3_4 = ay3_end * dt

    # Update positions and velocities
    x1 += (kx1_1 + 2 * kx1_2 + 2 * kx1_3 + kx1_4) / 6
    y1 += (ky1_1 + 2 * ky1_2 + 2 * ky1_3 + ky1_4) / 6
    vx1 += (kvx1_1 + 2 * kvx1_2 + 2 * kvx1_3 + kvx1_4) / 6
    vy1 += (kvy1_1 + 2 * kvy1_2 + 2 * kvy1_3 + kvy1_4) / 6

    x2 += (kx2_1 + 2 * kx2_2 + 2 * kx2_3 + kx2_4) / 6
    y2 += (ky2_1 + 2 * ky2_2 + 2 * ky2_3 + ky2_4) / 6
    vx2 += (kvx2_1 + 2 * kvx2_2 + 2 * kvx2_3 + kvx2_4) / 6
    vy2 += (kvy2_1 + 2 * kvy2_2 + 2 * kvy2_3 + kvy2_4) / 6

    x3 += (kx3_1 + 2 * kx3_2 + 2 * kx3_3 + kx3_4) / 6
    y3 += (ky3_1 + 2 * ky3_2 + 2 * ky3_3 + ky3_4) / 6
    vx3 += (kvx3_1 + 2 * kvx3_2 + 2 * kvx3_3 + kvx3_4) / 6
    vy3 += (kvy3_1 + 2 * kvy3_2 + 2 * kvy3_3 + kvy3_4) / 6

    return x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3

perturbation_factor = 1.5

# Assuming stars are on the x-axis initially and the planet is further out on the x-axis
r = 1.5 * pow(10, 11)  # Initial separation between stars, 1 AU for simplicity
r_planet = 3 * pow(10, 11) * perturbation_factor  # Initial distance of the planet from the system's center of mass

# Initial positions
x1, y1 = -r/2, 0  # Star 1
x2, y2 = r/2, 0   # Star 2
x3, y3 = r_planet, 0  # Planet

# Initial velocities for circular orbits
v_orbital = np.sqrt(G * (m1 + m2) / r)
v_planet = np.sqrt(G * (m1 + m2) / r_planet)

vx1, vy1 = 0, v_orbital * m2 / (m1 + m2)
vx2, vy2 = 0, -v_orbital * m1 / (m1 + m2)
vx3, vy3 = 0, v_planet

# Time setup
dt = 150  # Time step in seconds
total_time = 365.25 * 24 * 3600 * 10
times = np.arange(0, total_time, dt)

# Initialize arrays for plotting
x1s, y1s, x2s, y2s, x3s, y3s = [0], [0], [0], [0], [0], [0]
vx3s, vy3s = [], []

def semi_major_axis(x_initial, y_initial, x_final, y_final):
    # Calculate the distance between the initial and final positions
    distance_initial_final = np.sqrt((x_final - x_initial)**2 + (y_final - y_initial)**2)
    # Semi-major axis is half of the distance between the initial and final positions
    semi_major = distance_initial_final / 2
    # Convert semi-major axis from meters to Astronomical Units (AU)
    semi_major_AU = semi_major / (1.496 * 10**11)
    return semi_major_AU

# Call the function for each celestial body's orbit
semi_major_star1 = semi_major_axis(x1s[0], y1s[0], x1s[-1], y1s[-1])
semi_major_star2 = semi_major_axis(x2s[0], y2s[0], x2s[-1], y2s[-1])
semi_major_planet = semi_major_axis(x3s[0], y3s[0], x3s[-1], y3s[-1])


def orbital_period_years(semi_major_axis, m1, m2):
    G = 6.6726e-11  # m^3 kg^-1 s^-2
    seconds_in_year = 365.25 * 24 * 3600  # Approximation
    # the orbital period using Kepler's third law
    orbital_period_squared = (4 * np.pi**2 / (G * (m1 + m2))) * semi_major_axis**3
    # Take the square root to find the orbital period in seconds
    orbital_period_seconds = np.sqrt(orbital_period_squared)
    # Convert orbital period from seconds to years
    orbital_period_years = orbital_period_seconds / seconds_in_year
    return orbital_period_years

# Call the function for each celestial body's orbit
orbital_period_star1_years = orbital_period_years(semi_major_star1 * (1.496 * 10**11), m1, m2)  # Convert semi-major axis back to meters
orbital_period_star2_years = orbital_period_years(semi_major_star2 * (1.496 * 10**11), m1, m2)  # Convert semi-major axis back to meters
orbital_period_planet_years = orbital_period_years(semi_major_planet * (1.496 * 10**11), m1, m2)  # Convert semi-major axis back to meters


# Simulation loop
for _ in times:
    x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3 = runge_kutta(
        m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt
    )
    x1s.append(x1)
    y1s.append(y1)
    x2s.append(x2)
    y2s.append(y2)
    x3s.append(x3)
    y3s.append(y3)
    vx3s.append(vx3)
    vy3s.append(vy3)

# Calculate the separation of the planet from the center of mass
separation_center_of_mass = np.sqrt(np.array(x3s)**2 + np.array(y3s)**2)
velocity_integral_planet = [0]


# Plot velocity against separation
plt.figure(figsize=(10, 6))
plt.plot(separation_center_of_mass, velocity_integral_planet)
plt.xlabel('Separation from Center of Mass (m)')
plt.ylabel('Velocity of the Planet (m/s)')
plt.title('Velocity of the Planet vs. Separation from Center of Mass')

plt.tick_params(axis='both', direction='in', which='both', width=1, length=4, bottom=True, top=True, left=True, right=True)
plt.minorticks_on()

plt.show()
