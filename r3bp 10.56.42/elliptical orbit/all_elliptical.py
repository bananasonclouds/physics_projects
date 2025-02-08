import numpy as np
import matplotlib.pyplot as plt
import random

# Constants
G = 6.6726 * pow(10, -11)  # m^3 kg^-1 s^-2
M_sun= 1.99  * pow(10, 30) #in kg

# Masses of the two stars and the planet (in kg)
m1 = 4.95 * M_sun  # Star 1
m2 = 2.20*M_sun # Star 2
m3 = 2 * pow(10, 27)  #jupiter mass

miu= m2 / ( m1 + m2 ) #mass ratio
print(miu)


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

# Runge-Kutta function including the planet
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

    return x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3,ax1,ay1,ax2,ay2,ax3,ay3

# Run simulation including the planet
def run_runge(m1, m2, m3, t_max, dt):
    r = 1.5e11  # Initial distance between the stars (in meters) #2AU
    r_planet = 18e11 #20.589835e11 # * perturbation_factor_planet  # Initial distance of the planet from the center of mass OF BOTH STARS (in meters)
    
    # Velocities for stable orbits
    v_stars = np.sqrt(G * (m1 + m2) / r) * perturbation_factor
    v_planet = np.sqrt(G * (m1 + m2) / r_planet)  

    # Initial positions
    # Star 1
    x1, y1 = -m2 / (m1 + m2) * r, 0
    vx1, vy1 = 0, m2 / (m1 + m2) * v_stars
    # Star 2
    x2, y2 = m1 / (m1 + m2) * r, 0
    vx2, vy2 = 0, -m1 / (m1 + m2) * v_stars
    # Planet
    x3, y3 = 0, r_planet
    vx3, vy3 = -v_planet, 0

    x1s, y1s, x2s, y2s, x3s, y3s,vx1s, vy1s,vx2s,vy2s,vx3s,vy3s,ax1s,ay1s,ax2s,ay2s,ax3s,ay3s = [], [], [], [], [], [],[], [], [], [], [], [],[], [], [], [], [], []

    for _ in range(int(t_max / dt)):
        x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3,ax1,ay1,ax2,ay2,ax3,ay3 = runge_kutta(
            m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt
        )
        x1s.append(x1)
        y1s.append(y1)
        x2s.append(x2)
        y2s.append(y2)
        x3s.append(x3)
        y3s.append(y3)
        vx1s.append(vx1) 
        vy1s.append(vy1) 
        vx2s.append(vx2) 
        vy2s.append(vy2) 
        vx3s.append(vx3) 
        vy3s.append(vy3)
        ax1s.append(ax1) 
        ay1s.append(ay1) 
        ax2s.append(ax2) 
        ay2s.append(ay2) 
        ax3s.append(ax3) 
        ay3s.append(ay3)
    vx1s=np.array(vx1s)
    vy1s=np.array(vy1s)
    vx2s=np.array(vx2s)
    vy2s=np.array(vy2s)
    vx3s=np.array(vx3s)
    vy3s=np.array(vy3s)
    ax1s=np.array(ax1s)
    ay1s=np.array(ay1s)
    ax2s=np.array(ax2s)
    ay2s=np.array(ay2s)
    ax3s=np.array(ax3s)
    ay3s=np.array(ay3s)
    
    return  x1s, y1s, x2s, y2s, x3s, y3s,vx1s, vy1s,vx2s,vy2s,vx3s,vy3s,ax1s,ay1s,ax2s,ay2s,ax3s,ay3s

# Parameters for the simulation
dt = 250  # Time step
t_max = 365 * 24 * 3600 * 50
perturbation_factor = 0.5  
perturbation_factor_planet = 24 #5.59478  #generate_random_float() * 1
times=np.arange(0,t_max,dt)

# Run the simulation
x1s, y1s, x2s, y2s, x3s, y3s,vx1s, vy1s,vx2s,vy2s,vx3s,vy3s,ax1s,ay1s,ax2s,ay2s,ax3s,ay3s = run_runge(m1, m2, m3,t_max, dt)

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
orbital_period_star1_years = orbital_period_years(semi_major_star1 , m1, m2)  # Convert semi-major axis back to meters
orbital_period_star2_years = orbital_period_years(semi_major_star2 , m1, m2)  # Convert semi-major axis back to meters
orbital_period_planet_years = orbital_period_years(semi_major_planet , m1, m2)  # Convert semi-major axis back to meters

def calculate_separation(x1s, y1s, x2s, y2s, x3s, y3s):
    separation_star1 = np.sqrt((np.array(x3s) - np.array(x1s))**2 + (np.array(y3s) - np.array(y1s))**2)
    separation_star2 = np.sqrt((np.array(x3s) - np.array(x2s))**2 + (np.array(y3s) - np.array(y2s))**2)
    separation_COM = np.sqrt(np.array(x3s)**2 + np.array(y3s)**2)
    return separation_star1, separation_star2, separation_COM

separation_star1, separation_star2,separation_COM = calculate_separation(x1s, y1s, x2s, y2s, x3s, y3s)

accel_1=np.sqrt(pow(ax1s,2)+pow(ay1s,2))
accel_2=np.sqrt(pow(ax2s,2)+pow(ay2s,2))
accel_3=np.sqrt(pow(ax3s,2)+pow(ay3s,2))

vel_1=np.sqrt(pow(vx1s,2)+pow(vy1s,2))
vel_2=np.sqrt(pow(vx2s,2)+pow(vy2s,2))
vel_3=np.sqrt(pow(vx3s,2)+pow(vy3s,2))

'''
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

#velocity
ax1.plot(times, vel_3, label='Velocity (Planet)', color='green')
#ax1.plot(times, vel_2,label='star2 velocity',color='blue')
#ax1.plot(times, vel_1,label='star1 velocity',color='cyan')
ax1.tick_params(axis='both', direction='in', which='both', width=1, length=4, bottom=True, top=True, left=True, right=True)
ax1.minorticks_on()
ax1.set_ylabel('Velocity (ms$^{-1}$)')
#ax1.legend()

#acceleration
ax2.plot(times, accel_3, color='red',label='Acceleration (Planet)')
#ax2.plot(times, accel_2, label='Star 2 Acceleration', color='blue')
#ax2.plot(times, accel_1, label='Star 1 Acceleration', color='cyan')
ax2.set_ylabel('Acceleration (ms$^{-2}$)')
#ax2.legend()
ax2.tick_params(axis='both', direction='in', which='both', width=1, length=4, bottom=True, top=True, left=True, right=True)
ax2.minorticks_on()

#separation from stars
ax3.plot(times, separation_star1, label='Separation from Star 1', color='cyan', linestyle='-')
ax3.plot(times, separation_star2, label='Separation from Star 2', color='blue', linestyle='-')
ax3.plot(times,separation_COM, label="Separation from COM", color='darkblue',linestyle='--')
ax3.set_ylabel('Separation (m)')
ax3.set_xlabel('Time/(s)')
ax3.legend()
ax3.tick_params(axis='both', direction='in', which='both', width=1, length=4, bottom=True, top=True, left=True, right=True)
ax3.minorticks_on()'''

# Print the results
print("Orbital period of Star 1:", orbital_period_star1_years, "years")
print("Orbital period of Star 2:", orbital_period_star2_years, "years")
print("Orbital period of the Planet:", orbital_period_planet_years, "years")

# Print the results
print("Semi-major axis of Star 1:", semi_major_star1, "m")
print("Semi-major axis of Star 2:", semi_major_star2, "m")
print("Semi-major axis of the Planet:", semi_major_planet, "m")

plt.tight_layout()
plt.figure(figsize=(8, 6))
plt.plot(x1s, y1s, label='Star 1',color='blue')
plt.plot(x2s, y2s, label='Star 2',color='cyan')
plt.plot(x3s, y3s, label='Planet',color='red',linewidth=0.8)
plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
#plt.title('Binary Star System with P-type Orbit')
plt.axis('equal')
plt.tick_params(axis='both', direction='in', which='both', width=1, length=4, bottom=True, top=True, left=True, right=True)
plt.minorticks_on()
plt.scatter([x1s[0]], [y1s[0]], color='blue', marker='o', label='Star 1 Initial Position')
plt.scatter([x2s[0]], [y2s[0]], color='cyan', marker='o', label='Star 2 Initial Position')
plt.scatter([x3s[0]], [y3s[0]], color='red', marker='o', label='Planet Initial Position')
plt.show()

