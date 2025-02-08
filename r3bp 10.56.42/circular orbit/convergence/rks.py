import numpy as np
import matplotlib.pyplot as plt

G = 6.6726 * pow(10, -11)  # m^3 kg^-1 s^-2
M_sun= 1.99  * pow(10, 30) #in kg
M_jupiter= 2 * pow(10, 27) #jupiter mass

# Masses of the two stars and the planet (in kg)
m1 = 6 * M_sun  # Star 1
m2 = 2 * M_sun  # Star 2
m3 = 1

miu = m2 / (m1 + m2) # mass ratio

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

def runge_kutta2(m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt):
    ax1, ay1, ax2, ay2, ax3, ay3 = acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3)

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

    x1_mid = x1 + kx1_1 / 2
    y1_mid = y1 + ky1_1 / 2
    vx1_mid = vx1 + kvx1_1 / 2
    vy1_mid = vy1 + kvy1_1 / 2

    x2_mid = x2 + kx2_1 / 2
    y2_mid = y2 + ky2_1 / 2
    vx2_mid = vx2 + kvx2_1 / 2
    vy2_mid = vy2 + kvy2_1 / 2

    x3_mid = x3 + kx3_1 / 2
    y3_mid = y3 + ky3_1 / 2
    vx3_mid = vx3 + kvx3_1 / 2
    vy3_mid = vy3 + kvy3_1 / 2

    kx1_2 = vx1_mid * dt
    ky1_2 = vy1_mid * dt
    kvx1_2, kvy1_2, kvx2_2, kvy2_2, kvx3_2,kvy3_2  = acceleration(m1, m2, m3, x1_mid, y1_mid, x2_mid, y2_mid, x3_mid, y3_mid)

    kx2_2 = vx2_mid * dt
    ky2_2 = vy2_mid * dt

    kx3_2 = vx3_mid * dt
    ky3_2 = vy3_mid * dt

    x1 += kx1_2
    y1 += ky1_2
    vx1 += kvx1_2
    vy1 += kvy1_2

    x2 += kx2_2
    y2 += ky2_2
    # Note: We need to compute kvx2_2 and kvy2_2 here too if needed in further computations
    vx2 += kvx2_2
    vy2 += kvy2_2

    x3 += kx3_2
    y3 += ky3_2
    # Note: We need to compute kvx3_2 and kvy3_2 here too if needed in further computations
    vx3 += kvx3_2
    vy3 += kvy3_2

    return x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3

def runge_kutta3(m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt):
    ax1, ay1, ax2, ay2, ax3, ay3 = acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3)

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

    x1_mid = x1 + kx1_1 / 2
    y1_mid = y1 + ky1_1 / 2
    vx1_mid = vx1 + kvx1_1 / 2
    vy1_mid = vy1 + kvy1_1 / 2

    x2_mid = x2 + kx2_1 / 2
    y2_mid = y2 + ky2_1 / 2
    vx2_mid = vx2 + kvx2_1 / 2
    vy2_mid = vy2 + kvy2_1 / 2

    x3_mid = x3 + kx3_1 / 2
    y3_mid = y3 + ky3_1 / 2
    vx3_mid = vx3 + kvx3_1 / 2
    vy3_mid = vy3 + kvy3_1 / 2

    kx1_2 = vx1_mid * dt
    ky1_2 = vy1_mid * dt
    kvx1_2, kvy1_2, kvx2_2, kvy2_2, kvx3_2, kvy3_2 = acceleration(m1, m2, m3, x1_mid, y1_mid, x2_mid, y2_mid, x3_mid, y3_mid)

    x1_mid = x1 + kx1_2 / 2
    y1_mid = y1 + ky1_2 / 2
    vx1_mid = vx1 + kvx1_2 / 2
    vy1_mid = vy1 + kvy1_2 / 2

    kx2_2 = vx2_mid * dt
    ky2_2 = vy2_mid * dt

    x2_mid = x2 + kx2_2 / 2
    y2_mid = y2 + ky2_2 / 2
    vx2_mid = vx2 + kvx2_2 / 2
    vy2_mid = vy2 + kvy2_2 / 2


    kx1_3 = vx1_mid * dt
    ky1_3 = vy1_mid * dt
    kvx1_3, kvy1_3, kvx2_3, kvy2_3, kvx3_3, kvy3_3 = acceleration(m1, m2, m3, x1_mid, y1_mid, x2_mid, y2_mid, x3_mid, y3_mid)

    x1 += kx1_3
    y1 += ky1_3
    vx1 += kvx1_3
    vy1 += kvy1_3

    kx3_2 = vx3_mid * dt
    ky3_2 = vy3_mid * dt

    kx2_3 = vx3_mid * dt
    ky2_3 = vy3_mid * dt
    
    kx3_3 = (vx3 + kvx3_2 / 2) * dt
    ky3_3 = (vy3 + kvy3_2 / 2) * dt
    
    x3_mid = x3 + kx3_2 / 2
    y3_mid = y3 + ky3_2 / 2
    vx3_mid = vx3 + kvx3_2 / 2
    vy3_mid = vy3 + kvy3_2 / 2

    x2 += kx2_3
    y2 += ky2_3
    vx2 += kvx2_3
    vy2 += kvy2_3

    x3 += kx3_3
    y3 += ky3_3
    vx3 += kvx3_3
    vy3 += kvy3_3

    return x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3

def runge_kutta4(m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt):
    ax1, ay1, ax2, ay2, ax3, ay3 = acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3)

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

    x1_mid = x1 + kx1_1 / 2
    y1_mid = y1 + ky1_1 / 2
    x2_mid = x2 + kx2_1 / 2
    y2_mid = y2 + ky2_1 / 2
    x3_mid = x3 + kx3_1 / 2
    y3_mid = y3 + ky3_1 / 2

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

    x1_mid = x1 + kx1_2 / 2
    y1_mid = y1 + ky1_2 / 2
    x2_mid = x2 + kx2_2 / 2
    y2_mid = y2 + ky2_2 / 2
    x3_mid = x3 + kx3_2 / 2
    y3_mid = y3 + ky3_2 / 2

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

    x1_end = x1 + kx1_3
    y1_end = y1 + ky1_3
    vx1_end = vx1 + kvx1_3
    vy1_end = vy1 + kvy1_3
    x2_end = x2 + kx2_3
    y2_end = y2 + ky2_3
    vx2_end = vx2 + kvx2_3
    vy2_end = vy2 + kvy2_3
    x3_end = x3 + kx3_3
    y3_end = y3 + ky3_3
    vx3_end = vx3 + kvx3_3
    vy3_end = vy3 + kvy3_3

    return x1_end, y1_end, vx1_end, vy1_end, x2_end, y2_end, vx2_end, vy2_end, x3_end, y3_end, vx3_end, vy3_end
def run_runge(m1, m2, m3, times, dt, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3):
    x1s_rk2, y1s_rk2, x2s_rk2, y2s_rk2, x3s_rk2, y3s_rk2 = [], [], [], [], [], []
    x1s_rk3, y1s_rk3, x2s_rk3, y2s_rk3, x3s_rk3, y3s_rk3 = [], [], [], [], [], []
    x1s_rk4, y1s_rk4, x2s_rk4, y2s_rk4, x3s_rk4, y3s_rk4 = [], [], [], [], [], []

    for _ in times:
        x1_rk2, y1_rk2, vx1_rk2, vy1_rk2, x2_rk2, y2_rk2, vx2_rk2, vy2_rk2, x3_rk2, y3_rk2, vx3_rk2, vy3_rk2 = runge_kutta2(
            m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt
        )
        x1s_rk2.append(x1_rk2)
        y1s_rk2.append(y1_rk2)
        x2s_rk2.append(x2_rk2)
        y2s_rk2.append(y2_rk2)
        x3s_rk2.append(x3_rk2)
        y3s_rk2.append(y3_rk2)
        
        x1_rk3, y1_rk3, vx1_rk3, vy1_rk3, x2_rk3, y2_rk3, vx2_rk3, vy2_rk3, x3_rk3, y3_rk3, vx3_rk3, vy3_rk3 = runge_kutta3(
            m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt
        )
        x1s_rk3.append(x1_rk3)
        y1s_rk3.append(y1_rk3)
        x2s_rk3.append(x2_rk3)
        y2s_rk3.append(y2_rk3)
        x3s_rk3.append(x3_rk3)
        y3s_rk3.append(y3_rk3)

        x1_rk4, y1_rk4, vx1_rk4, vy1_rk4, x2_rk4, y2_rk4, vx2_rk4, vy2_rk4, x3_rk4, y3_rk4, vx3_rk4, vy3_rk4 = runge_kutta3(
            m1, m2, m3, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3, dt)
        x1s_rk4.append(x1_rk4)
        y1s_rk4.append(y1_rk4)
        x2s_rk4.append(x2_rk4)
        y2s_rk4.append(y2_rk4)
        x3s_rk4.append(x3_rk4)
        y3s_rk4.append(y3_rk4)

    x1s_rk2 = np.array(x1s_rk2)
    y1s_rk2 = np.array(y1s_rk2)
    x2s_rk2 = np.array(x2s_rk2)
    y2s_rk2 = np.array(y2s_rk2)
    x3s_rk2 = np.array(x3s_rk2)
    y3s_rk2 = np.array(y3s_rk2)
    
    x1s_rk3 = np.array(x1s_rk3)
    y1s_rk3 = np.array(y1s_rk3)
    x2s_rk3 = np.array(x2s_rk3)
    y2s_rk3 = np.array(y2s_rk3)
    x3s_rk3 = np.array(x3s_rk3)
    y3s_rk3 = np.array(y3s_rk3)

    x1s_rk4 = np.array(x1s_rk4)
    y1s_rk4 = np.array(y1s_rk4)
    x2s_rk4 = np.array(x2s_rk4)
    y2s_rk4 = np.array(y2s_rk4)
    x3s_rk4 = np.array(x3s_rk4)
    y3s_rk4 = np.array(y3s_rk4)

    return (
        x1s_rk2, y1s_rk2, x2s_rk2, y2s_rk2, x3s_rk2, y3s_rk2,
        x1s_rk3, y1s_rk3, x2s_rk3, y2s_rk3, x3s_rk3, y3s_rk3,
        x1s_rk4, y1s_rk4, x2s_rk4, y2s_rk4, x3s_rk4, y3s_rk4
    )

# Run the simulation
t_max = 365.25 * 24 * 3600 * 3
dt = 200
times = np.arange(0, t_max, dt)
r = 1.5 * pow(10, 11)  # Initial separation between stars, 1 AU for simplicity
r_planet = 7.5 * pow(10, 11) #9.3AU  # Initial distance of the planet from the system's center of mass
v_orbital = np.sqrt(G * (m1 + m2) / r) 
v_planet = np.sqrt(G * (m1 + m2) / r_planet)

# Define initial conditions
x1_initial, y1_initial = -r/2, 0
x2_initial, y2_initial = r/2, 0
x3_initial, y3_initial = r_planet, 0
vx1_initial, vy1_initial = 0, v_orbital* m2 / (m1 + m2)
vx2_initial, vy2_initial = 0,-v_orbital * m1 / (m1 + m2) 
vx3_initial, vy3_initial = 0, v_planet

(
    x1s_rk2, y1s_rk2, x2s_rk2, y2s_rk2, x3s_rk2, y3s_rk2,
    x1s_rk3, y1s_rk3, x2s_rk3, y2s_rk3, x3s_rk3, y3s_rk3,
    x1s_rk4, y1s_rk4, x2s_rk4, y2s_rk4, x3s_rk4, y3s_rk4
) = run_runge(
    m1, m2, m3, times, dt,
    x1_initial, y1_initial, vx1_initial, vy1_initial,
    x2_initial, y2_initial, vx2_initial, vy2_initial,
    x3_initial, y3_initial, vx3_initial, vy3_initial,
)

# Plotting
plt.figure(figsize=(12, 8))
plt.title('Orbits of Stars and Planet')

plt.plot(x1s_rk2, y1s_rk2, label='Star 1 (RK2)', color='blue')
plt.plot(x2s_rk2, y2s_rk2, label='Star 2 (RK2)', color='red')
plt.plot(x3s_rk2, y3s_rk2, label='Planet (RK2)', color='green')

plt.plot(x1s_rk3, y1s_rk3, label='Star 1 (RK3)', color='blue', linestyle='--')
plt.plot(x2s_rk3, y2s_rk3, label='Star 2 (RK3)', color='red', linestyle='--')
plt.plot(x3s_rk3, y3s_rk3, label='Planet (RK3)', color='green', linestyle='--')

plt.plot(x1s_rk4, y1s_rk4, label='Star 1 (RK4)', color='blue', linestyle='--')
plt.plot(x2s_rk4, y2s_rk4, label='Star 2 (RK4)', color='red', linestyle='--')
plt.plot(x3s_rk4, y3s_rk4, label='Planet (RK4)', color='green', linestyle='--')


plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.legend()
plt.grid(True)
plt.show()
