import numpy as np
import matplotlib.pyplot as plt

me = 5.9742 * pow(10, 24)  # Mass of the Earth in kg
mm = 7.35 * pow(10, 22)  # Mass of the Moon in kg
G = 6.6726 * pow(10, -11)  # Gravitational constant in m^3 kg^-1 s^-2
r12 = 384400.0 * 10 ** 3  # Distance between the CoM of the Earth and Moon in meters
M = me + mm

re = r12 * mm / (mm + me)
rm = r12 * me / (mm + me)

m_tot = me + mm
delta_r = r12 * pow(mm / (3 * me), 1 / 3)
T = 27.3*24*3600  #np.sqrt(4 * pow(np.pi, 2) * pow(r12, 3) / (G * m_tot))

ts = np.linspace(0, 2400000, 100000)
rocket_xs = []
rocket_ys = []
l2_xs = []
l2_ys = []
distance_values = []  # List to store distance values

for i in range(len(ts)):
    t = ts[i]
    rocket_x = r12 * np.cos(2 * np.pi * t / T) + delta_r - re
    rocket_y = r12 * np.sin(2 * np.pi * t / T)
    l2_x = -re * np.cos(2 * np.pi * t / T)
    l2_y = -re * np.sin(2 * np.pi * t / T)
    
    rocket_xs.append(rocket_x)
    rocket_ys.append(rocket_y)
    l2_xs.append(l2_x)
    l2_ys.append(l2_y)
    
    # Calculate the distance between rocket and L2 point
    distance = np.sqrt((rocket_x - l2_x) ** 2 + (rocket_y - l2_y) ** 2) -r12
    distance_values.append(distance)

plt.figure(figsize=(10, 5))
plt.plot(ts, distance_values)
plt.title("Distance between Rocket and L2 Point")
plt.xlabel("Time (seconds)")
plt.ylabel("Distance (meters)")

plt.tight_layout()
plt.show()
