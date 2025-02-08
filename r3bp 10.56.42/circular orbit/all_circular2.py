import numpy as np
import matplotlib.pyplot as plt
import math

G = 6.6726 * pow(10, -11)  # m^3 kg^-1 s^-2
M_sun= 1.99  * pow(10, 30) #in kg

# Masses of the two stars and the planet (in kg)
m1 = 6 * M_sun  # Star 1
m2 = 6 *M_sun # Star 2
m3 = 1

miu= m2 / ( m1 + m2 ) #mass ratio
print(miu)

# Acceleration function including the planet
def acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3):
    r12 = ((x2 - x1)**2 + (y2 - y1)**2)

    theta=math.atan2(y2-y1, x2-x1)
    x1_acc=G*m2/r12*np.cos(theta)
    y1_acc=G*m2/r12*np.sin(theta)
    x2_acc=-G*m1/r12*np.cos(theta)
    y2_acc=-G*m1/r12*np.sin(theta)

    x3_acc= G * m1 * (x1-x3)/(pow((pow(x3-x1,2)+pow(y3-y1,2)),3/2))  +G *m2*(x2-x3)/(pow(pow(x2-x3,2)+pow(y2-y3,2),3/2))
    y3_acc= G * m1 * (y1-y3)/(pow((pow(x3-x1,2)+pow(y3-y1,2)),3/2))  +G *m2*(y2-y3)/(pow(pow(x2-x3,2)+pow(y2-y3,2),3/2))
    return x1_acc, y1_acc, x2_acc, y2_acc, x3_acc, y3_acc

def runge_kutta(x1,y1,x2,y2,x3,y3,x1_vel,y1_vel,x2_vel,y2_vel,x3_vel,y3_vel,x1_acc,y1_acc,x2_acc,y2_acc,x3_acc,y3_acc,a):
    
    x1_acc, y1_acc, x2_acc, y2_acc, x3_acc, y3_acc = acceleration(m1, m2, m3, x1, y1, x2, y2, x3, y3)

    z1x_1 = x1 + a / 2 * x1_vel
    z1y_1 = y1 + a / 2 * y1_vel

    z1x_2 = x2 + a / 2 * x2_vel
    z1y_2 = y2 + a / 2 * y2_vel

    z1x_3 = x3 + a / 2 * x3_vel
    z1y_3 = y3 + a / 2 * y3_vel

    z1x1_vel=x1_vel+a/2*x1_acc
    z1y1_vel=y1_vel+a/2*y1_acc

    z1x2_vel=x2_vel+a/2*x2_acc
    z1y2_vel=y2_vel+a/2*y2_acc

    z1x3_vel=x3_vel+a/2*x3_acc
    z1y3_vel=y3_vel+a/2*y3_acc

    z1x1_acc,z1y1_acc=acceleration(z1x_1,z1y_1)
    z1x2_acc,z1y2_acc=acceleration(z1x_2,z1y_2)
    z1x3_acc,z1y3_acc=acceleration(z1x_3,z1y_3)

    z2x_1=x1+a/2*z1x1_vel
    z2y_1=y1+a/2*z1y1_vel

    z2x_2=x2+a/2*z1x2_vel
    z2y_2=y2+a/2*z1y2_vel

    z2x_3=x3+a/2*z1x3_vel
    z2y_3=y3+a/2*z1y3_vel

    z2x1_vel=x1_vel+a/2*z1x1_acc
    z2y1_vel=y1_vel+a/2*z1y1_acc

    z2x2_vel=x2_vel+a/2*z1x2_acc
    z2y2_vel=y2_vel+a/2*z1y2_acc

    z2x3_vel=x3_vel+a/2*z1x3_acc
    z2y3_vel=y3_vel+a/2*z1y3_acc

    z2x1_acc,z2y1_acc=acceleration(z2x_1,z2y_1)
    z2x2_acc,z2y2_acc=acceleration(z2x_2,z2y_2)
    z2x3_acc,z2y3_acc=acceleration(z2x_3,z2y_3)

    z3x_1=x1+a*z2x1_vel
    z3y_1=y1+a*z2y1_vel

    z3x_2=x2+a*z2x2_vel
    z3y_2=y2+a*z2y2_vel

    z3x_3=x3+a*z2x3_vel
    z3y_3=y3+a*z2y3_vel

    z3x1_vel=x1_vel+a*z2x1_acc
    z3y1_vel=y1_vel+a*z2y1_acc

    z3x2_vel=x2_vel+a*z2x2_acc
    z3y2_vel=y2_vel+a*z2y2_acc

    z3x3_vel=x3_vel+a*z2x3_acc
    z3y3_vel=y3_vel+a*z2y3_acc

    z3x1_acc,z3y1_acc=acceleration(z3x_1,z3y_1)
    z3x2_acc,z3y2_acc=acceleration(z3x_2,z3y_2)
    z3x3_acc,z3y3_acc=acceleration(z3x_3,z3y_3)

    x_1=x1+(a/6)*(x1_vel+2*z1x1_vel+2*z2x1_vel+z3x1_vel)
    y_1=y1+(a/6)*(y1_vel+2*z1y1_vel+2*z2y1_vel+z3y1_vel)

    x_2=x2+(a/6)*(x2_vel+2*z1x2_vel+2*z2x2_vel+z3x2_vel)
    y_2=y2+(a/6)*(y2_vel+2*z1y2_vel+2*z2y2_vel+z3y2_vel)

    x_3=x3+(a/6)*(x3_vel+2*z1x3_vel+2*z2x3_vel+z3x3_vel)
    y_3=y3+(a/6)*(y3_vel+2*z1y3_vel+2*z2y3_vel+z3y3_vel)

    x1_vel=x1_vel+(a/6)*(x1_acc+2*z1x1_acc+2*z2x1_acc+z3x1_acc)
    y1_vel=y1_vel+(a/6)*(y1_acc+2*z1y1_acc+2*z2y1_acc+z3y1_acc)

    x2_vel=x2_vel+(a/6)*(x2_acc+2*z1x2_acc+2*z2x2_acc+z3x2_acc)
    y2_vel=y2_vel+(a/6)*(y2_acc+2*z1y2_acc+2*z2y2_acc+z3y2_acc)

    x3_vel=x3_vel+(a/6)*(x3_acc+2*z1x3_acc+2*z2x3_acc+z3x3_acc)
    y3_vel=y3_vel+(a/6)*(y3_acc+2*z1y3_acc+2*z2y3_acc+z3y3_acc)

    x1_acc,y1_acc=acceleration(x_1,y1)
    x2_acc,y2_acc=acceleration(x_2,y1)
    x3_acc,y3_acc=acceleration(x_3,y_3)

    return x_1,y_1,x1_vel,y1_vel,x1_acc,y1_acc  ,x_2,y_2,x2_vel,y2_vel,x2_acc,y2_acc  ,x_3,y_3,x3_vel,y3_vel,x3_acc,y3_acc

def run_runge(m1,m2,m3,times,dt):
    r = 1.5 * pow(10, 11)  # Initial separation between stars, 1 AU for simplicity
    r_planet = 14 * pow(10, 11) #9.3AU  # Initial distance of the planet from the system's center of mass

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
        x1, y1, vx1, vy1,ax1,ay1 = runge_kutta(
            m1, m2, m3, x1, y1, vx1, vy1, dt
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
    
    return  x1s, y1s, x2s, y2s, x3s, y3s,vx1s, vy1s,vx2s,vy2s,vx3s,vy3s,ax1s,ay1s,ax2s,ay2s,ax3s,ay3s,r,r_planet
 
dt = 50  # Time step in seconds
t_max = 365.25 * 24 * 3600 * 8 # One year in seconds
times = np.arange(0, t_max, dt)

x1s, y1s, x2s, y2s, x3s, y3s,vx1s, vy1s,vx2s,vy2s,vx3s,vy3s,ax1s,ay1s,ax2s,ay2s,ax3s,ay3s,r,r_planet = run_runge(m1, m2, m3,times, dt)

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
