import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.6726 * pow(10, -11)  #m^3 kg^-1 s^-2

# Masses of the two stars (in kg)
m1 = 0.45 * 1.99*pow(10,30) #1.5 * 1.989 * pow(10, 30)  
m2 = 0.2*1.99 * pow(10, 30) 

def acceleration(m1, m2, x1, y1, x2, y2):
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    F = G * m1 * m2 / r**2
    Fx = F * (x2 - x1) / r
    Fy = F * (y2 - y1) / r
    return Fx / m1, Fy / m1, -Fx / m2, -Fy / m2

def runge_kutta(m1, m2, x1, y1, vx1, vy1, x2, y2, vx2, vy2, dt):
    # First set 
    ax1, ay1, ax2, ay2 = acceleration(m1, m2, x1, y1, x2, y2)
    
    kx1_1 = vx1 * dt
    ky1_1 = vy1 * dt
    kvx1_1 = ax1 * dt
    kvy1_1 = ay1 * dt

    kx2_1 = vx2 * dt
    ky2_1 = vy2 * dt
    kvx2_1 = ax2 * dt
    kvy2_1 = ay2 * dt
    
    # Second set 
    ax1_mid, ay1_mid, ax2_mid, ay2_mid = acceleration(
        m1, m2, 
        x1 + kx1_1 / 2, y1 + ky1_1 / 2, 
        x2 + kx2_1 / 2, y2 + ky2_1 / 2
    )
    
    kx1_2 = (vx1 + kvx1_1 / 2) * dt
    ky1_2 = (vy1 + kvy1_1 / 2) * dt
    kvx1_2 = ax1_mid * dt
    kvy1_2 = ay1_mid * dt

    kx2_2 = (vx2 + kvx2_1 / 2) * dt
    ky2_2 = (vy2 + kvy2_1 / 2) * dt
    kvx2_2 = ax2_mid * dt
    kvy2_2 = ay2_mid * dt

    # Third set 
    ax1_mid, ay1_mid, ax2_mid, ay2_mid = acceleration(
        m1, m2, 
        x1 + kx1_2 / 2, y1 + ky1_2 / 2, 
        x2 + kx2_2 / 2, y2 + ky2_2 / 2
    )
    
    kx1_3 = (vx1 + kvx1_2 / 2) * dt
    ky1_3 = (vy1 + kvy1_2 / 2) * dt
    kvx1_3 = ax1_mid * dt
    kvy1_3 = ay1_mid * dt

    kx2_3 = (vx2 + kvx2_2 / 2) * dt
    ky2_3 = (vy2 + kvy2_2 / 2) * dt
    kvx2_3 = ax2_mid * dt
    kvy2_3 = ay2_mid * dt

    # Fourth set 
    ax1_end, ay1_end, ax2_end, ay2_end = acceleration(
        m1, m2, 
        x1 + kx1_3, y1 + ky1_3, 
        x2 + kx2_3, y2 + ky2_3
    )
    
    kx1_4 = (vx1 + kvx1_3) * dt
    ky1_4 = (vy1 + kvy1_3) * dt
    kvx1_4 = ax1_end * dt
    kvy1_4 = ay1_end * dt

    kx2_4 = (vx2 + kvx2_3) * dt
    ky2_4 = (vy2 + kvy2_3) * dt
    kvx2_4 = ax2_end * dt
    kvy2_4 = ay2_end * dt

    x1 += (kx1_1 + 2 * kx1_2 + 2 * kx1_3 + kx1_4) / 6
    y1 += (ky1_1 + 2 * ky1_2 + 2 * ky1_3 + ky1_4) / 6
    vx1 += (kvx1_1 + 2 * kvx1_2 + 2 * kvx1_3 + kvx1_4) / 6
    vy1 += (kvy1_1 + 2 * kvy1_2 + 2 * kvy1_3 + kvy1_4) / 6

    x2 += (kx2_1 + 2 * kx2_2 + 2 * kx2_3 + kx2_4) / 6
    y2 += (ky2_1 + 2 * ky2_2 + 2 * ky2_3 + ky2_4) / 6
    vx2 += (kvx2_1 + 2 * kvx2_2 + 2 * kvx2_3 + kvx2_4) / 6
    vy2 += (kvy2_1 + 2 * kvy2_2 + 2 * kvy2_3 + kvy2_4) / 6

    return x1, y1, vx1, vy1, x2, y2, vx2, vy2

def run_runge(m1, m2, perturbation_factor, t_max, dt):
    r = 1.5e11  #meters
    v = np.sqrt(G * (m1 + m2) / r) * perturbation_factor

    # Initial conditions centered around the center of mass
    x1, y1 = -m2 / (m1 + m2) * r, 0
    x2, y2 = m1 / (m1 + m2) * r, 0
    vx1, vy1 = 0, m2 / (m1 + m2) * v
    vx2, vy2 = 0, -m1 / (m1 + m2) * v

    x1s, y1s, x2s, y2s = [], [], [], []

    for _ in range(int(t_max / dt)):
        x1, y1, vx1, vy1, x2, y2, vx2, vy2 = runge_kutta(m1, m2, x1, y1, vx1, vy1, x2, y2, vx2, vy2, dt)
        x1s.append(x1)
        y1s.append(y1)
        x2s.append(x2)
        y2s.append(y2)

    return x1s, y1s, x2s, y2s

#parameters
t_max = 365 * 24 * 3600  
dt = 200
perturbation_factor = 0.5 
x1s, y1s, x2s, y2s = run_runge(m1, m2, perturbation_factor, t_max, dt)

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(x1s, y1s, label='Star 1',color='blue')
plt.plot(x2s, y2s, label='Star 2',color='cyan')
plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.legend()
#plt.title('Binary Star System ')
plt.axis('equal')  
plt.show()
