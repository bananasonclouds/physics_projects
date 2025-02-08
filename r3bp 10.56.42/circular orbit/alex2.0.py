import numpy as np
import matplotlib.pyplot as plt

G = 6.6726 * pow(10, -11)  # m^3 kg^-1 s^-2
M_sun= 1.99  * pow(10, 30) #in kg
M_jupiter= 2 * pow(10, 27) #jupiter mass

# Masses of the two stars and the planet (in kg)
m1 = 6 * M_sun  # Star 1
m2 = 2 *M_sun # Star 2
m3 = 0.01 * M_jupiter

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
def runge_kutta(x,y,x_vel,y_vel,x_acc,y_acc,a,t):
    z1x = x + a / 2 * x_vel
    z1y = y + a / 2 * y_vel
    z1x_vel=x_vel+a/2*x_acc
    z1y_vel=y_vel+a/2*y_acc
    z1x_acc,z1y_acc=acceleration(z1x,z1y,t+a/2)
    z2x=x+a/2*z1x_vel
    z2y=y+a/2*z1y_vel
    z2x_vel=x_vel+a/2*z1x_acc
    z2y_vel=y_vel+a/2*z1y_acc
    z2x_acc,z2y_acc=acceleration(z2x,z2y,t + a/2)
    z3x=x+a*z2x_vel
    z3y=y+a*z2y_vel
    z3x_vel=x_vel+a*z2x_acc
    z3y_vel=y_vel+a*z2y_acc
    z3x_acc,z3y_acc=acceleration(z3x,z3y,t+a)
    x1=x+(a/6)*(x_vel+2*z1x_vel+2*z2x_vel+z3x_vel)
    y1=y+(a/6)*(y_vel+2*z1y_vel+2*z2y_vel+z3y_vel)
    x1_vel=x_vel+(a/6)*(x_acc+2*z1x_acc+2*z2x_acc+z3x_acc)
    y1_vel=y_vel+(a/6)*(y_acc+2*z1y_acc+2*z2y_acc+z3y_acc)
    x1_acc,y1_acc=acceleration(x1,y1,t+a)
    
    return x1,y1,x1_vel,y1_vel,x1_acc,y1_acc


def run_runge(m1, m2, m3,times, dt):
    r = 0.5 * pow(10, 8)  # Initial separation between stars, 1 AU for simplicity
    r_planet = 14 * pow(10, 11) # perturbation_factor  # Initial distance of the planet from the system's center of mass

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

    x1s, y1s, x2s, y2s, x3s, y3s,vx1s, vy1s,vx2s,vy2s,vx3s,vy3s,ax1s,ay1s,ax2s,ay2s,ax3s,ay3s = [], [], [], [], [], [],[], [], [], [], [], [],[], [], [], [], [], []
    
    for _ in times:
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

# Time setup
dt = 250  # Time step in seconds
t_max = 365.25 * 24 * 3600 * 50 # One year in seconds
times = np.arange(0, t_max, dt)
perturbation_factor = 2 #1.5

x1s, y1s, x2s, y2s, x3s, y3s,vx1s, vy1s,vx2s,vy2s,vx3s,vy3s,ax1s,ay1s,ax2s,ay2s,ax3s,ay3s = run_runge(m1, m2, m3,times, dt)

def calculate_separation(x1s, y1s, x2s, y2s, x3s, y3s):
    separation_star1 = np.sqrt((np.array(x3s) - np.array(x1s))**2 + (np.array(y3s) - np.array(y1s))**2)
    separation_star2 = np.sqrt((np.array(x3s) - np.array(x2s))**2 + (np.array(y3s) - np.array(y2s))**2)
    separation_COM = np.sqrt(np.array(x3s)**2 + np.array(y3s)**2)
    return separation_star1, separation_star2 ,separation_COM

separation_star1, separation_star2,separation_COM = calculate_separation(x1s, y1s, x2s, y2s, x3s, y3s)

# Create subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

accel_1=np.sqrt(pow(ax1s,2)+pow(ay1s,2))
accel_2=np.sqrt(pow(ax2s,2)+pow(ay2s,2))
accel_3=np.sqrt(pow(ax3s,2)+pow(ay3s,2))

vel_1=np.sqrt(pow(vx1s,2)+pow(vy1s,2))
vel_2=np.sqrt(pow(vx2s,2)+pow(vy2s,2))
vel_3=np.sqrt(pow(vx3s,2)+pow(vy3s,2))

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
ax3.minorticks_on()

'''# separation of planet from COM against time
ax4.plot(times, separation_COM, label="Separation from COM", color='green')
ax4.set_ylabel("Separation from COM (m)")
ax4.set_xlabel('Time/(s)')
ax4.tick_params(axis='both', direction='in', which='both', width=1, length=4, bottom=True, top=True, left=True, right=True)
ax4.minorticks_on()
ax4.legend()'''

plt.figure(figsize=(10, 6))
plt.plot(x1s, y1s, label="Star 1",color='cyan')
plt.plot(x2s, y2s, label="Star 2",color='blue')
plt.plot(x3s, y3s, label="Planet",color='red')
plt.xlabel("X Position (m)")
plt.ylabel("Y Position (m)")
plt.legend()
#plt.title("Simulation of a Binary Star System with a Planet")
plt.axis('equal')
plt.tick_params(axis='both', direction='in', which='both', width=1, length=4, bottom=True, top=True, left=True, right=True)
plt.minorticks_on()
plt.scatter([x1s[0]], [y1s[0]], color='cyan', marker='o', label='Star 1 Initial Position')
plt.scatter([x2s[0]], [y2s[0]], color='blue', marker='o', label='Star 2 Initial Position')
plt.scatter([x3s[0]], [y3s[0]], color='red', marker='o', label='Planet Initial Position')
plt.show()

