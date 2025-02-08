import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.6726 * pow(10, -11)  # m^3 kg^-1 s^-2
M_sun= 1.99  * pow(10, 30) #in kg

# Masses of the two stars and the planet (in kg)
m1 = 2 * M_sun  # Star 1
m2 = M_sun # Star 2
m3 = 5 * pow(10, 24)  # Planet (Earth mass) 5x10^24


# Acceleration function including the planet
def acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3):
    r12 = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    r13 = np.sqrt((x3 - x1)**2 + (y3 - y1)**2)
    r23 = np.sqrt((x3 - x2)**2 + (y3 - y2)**2)

    F12 = G * m1 * m2 / r12**2
    F13 = G * m1 * m3 / r13**2
    F23 = G * m2 * m3 / r23**2

    Fx12 = F12 * (x2 - x1) / r12
    Fy12 = F12 * (y2 - y1) / r12
    Fx13 = F13 * (x3 - x1) / r13
    Fy13 = F13 * (y3 - y1) / r13
    Fx23 = F23 * (x3 - x2) / r23
    Fy23 = F23 * (y3 - y2) / r23

    ax1 = (Fx12 + Fx13) / m1
    ay1 = (Fy12 + Fy13) / m1
    ax2 = (-Fx12 + Fx23) / m2
    ay2 = (-Fy12 + Fy23) / m2
    ax3 = (-Fx13 - Fx23) / m3
    ay3 = (-Fy13 - Fy23) / m3

    return ax1, ay1, ax2, ay2, ax3, ay3


def runge_kutta(m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt):
    ax1, ay1, ax2, ay2, ax3, ay3 = acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3)

# First set
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

    # Second set
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

    # Third set
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

    # Fourth set
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

def run_runge(m1, m2, m3, perturbation_factor, t_max, dt):
    r = 1e11   # Initial distance between the stars (in meters) 1e11 
    r_planet = 1.5 * 10 ** 9   #2 * 10 ** 9 # Initial distance of the planet from Star 1 (in meters)
    
    # Velocities for stable orbits
    v_stars = np.sqrt(G * (m1 + m2) / r) * perturbation_factor
    v_planet = np.sqrt(G * m1 / r_planet) * perturbation_factor_planet #?

    # Initial conditions
    # Star 1
    x1, y1 = -m2 / (m1 + m2) * r, 0
    vx1, vy1 = 0, m2 / (m1 + m2) * v_stars
    # Star 2
    x2, y2 = m1 / (m1 + m2) * r, 0
    vx2, vy2 = 0, -m1 / (m1 + m2) * v_stars
    # Planet orbiting Star 1
    x3, y3 = x1, y1 + r_planet
    vx3, vy3 = vx1 - np.sqrt(G * m1 / r_planet), vy1

    x1s, y1s, x2s, y2s, x3s, y3s = [], [], [], [], [], []

    # Run the simulation
    for _ in range(int(t_max / dt)):
        x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3 = runge_kutta(
            m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt
        )
        x1s.append(x1)
        y1s.append(y1)
        x2s.append(x2)
        y2s.append(y2)
        x3s.append(x3)
        y3s.append(y3)

    return x1s, y1s, x2s, y2s, x3s, y3s

# Parameters for the simulation
t_max = 365 * 24 * 3600 * 1  # Ten years
dt = 250  # Time step, longest = 20 dt
perturbation_factor = 0.5 # initial velocity perturbation factor for stable orbits #ideal pert factor is 0.5
perturbation_factor_planet = 0 #2.35

# Run the simulation
x1s, y1s, x2s, y2s, x3s, y3s = run_runge(m1, m2, m3, perturbation_factor, t_max, dt)

# Find the limits for the plot
x_min = min(min(x1s), min(x2s)) - 0.5 * abs(min(min(x1s), min(x2s)))
x_max = max(max(x1s), max(x2s)) + 0.5 * abs(max(max(x1s), max(x2s)))
y_min = min(min(y1s), min(y2s)) - 0.5 * abs(min(min(y1s), min(y2s)))
y_max = max(max(y1s), max(y2s)) + 0.5 * abs(max(max(y1s), max(y2s)))

# Plotting with adjusted axis limits
plt.figure(figsize=(8, 6))
plt.plot(x1s, y1s, label='Star 1')
plt.plot(x2s, y2s, label='Star 2')
plt.plot(x3s, y3s, label='Planet', linestyle='dotted',linewidth=3.5)
plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
plt.title('Binary Star System with S-type Planet')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.axis('equal')
plt.show()

 #orbit from center of mass

'''def run_runge(m1, m2, m3, perturbation_factor, t_max, dt):
    r = 1e11  # Initial distance between the stars (in meters)
    r_planet = 2e10  # Distance of the planet from the center of mass (in meters)
    
    # Velocities for stable orbits
    v_stars = np.sqrt(G * (m1 + m2) / r) * perturbation_factor
    v_planet = np.sqrt(G * (m1 + m2) / r_planet) * perturbation_factor

    # Initial conditions
    # Star 1
    x1, y1 = -m2 / (m1 + m2) * r, 0
    vx1, vy1 = 0, m2 / (m1 + m2) * v_stars
    # Star 2
    x2, y2 = m1 / (m1 + m2) * r, 0
    vx2, vy2 = 0, -m1 / (m1 + m2) * v_stars
    # Center of mass of the binary system
    com_x = (m1 * x1 + m2 * x2) / (m1 + m2)
    com_y = (m1 * y1 + m2 * y2) / (m1 + m2)
    # Planet
    x3, y3 = com_x, com_y + r_planet  # Positioning the planet at r_planet distance from the COM
    vx3, vy3 = -v_planet, 0  # Initial velocity of the planet

    '''
