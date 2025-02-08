import numpy as np
import matplotlib.pyplot as plt

me = 5.9742 *pow(10 , 24)  # Mass of the Earth in kg
mm = 7.35 * pow(10,22)  # Mass of the Moon in kg
G = 6.6726 * pow(10, -11)  # Gravitational constant in m^3 kg^-1 s^-2
r12 = 384400.0 * 10 ** 3  # Distance between the CoM of the Earth and Moon in meters
M = me + mm

re=r12*mm/(mm+me)
rm=r12*me/(mm+me)
#print(re)
m_tot=me+mm
delta_r=r12*pow(mm/(3*me),1/3)
T=np.sqrt(4*pow(np.pi,2)*pow(r12,3)/(G*m_tot))

ts=np.linspace(0,2400000,100000)
moon_xs=[]
moon_ys=[]
for i in range(len(ts)):
    t=ts[i]
    moon_x=rm*np.cos(2*np.pi*t/T)
    moon_y=rm*np.sin(2*np.pi*t/T)
    moon_xs.append(moon_x)
    moon_ys.append(moon_y)
plt.plot(moon_xs,moon_ys,linestyle='dashed',color='black')

#print(T)

perturbation_factor=1
v_init=2*np.pi*(r12+delta_r-re)/T*perturbation_factor

#v_init=50
def acceleration(x,y,t):
    earth_x=-re*np.cos(2*np.pi*t/T)
    earth_y=-re*np.sin(2*np.pi*t/T)
    moon_x=rm*np.cos(2*np.pi*t/T)
    moon_y=rm*np.sin(2*np.pi*t/T)
    x_acc= -G * me * (x-earth_x)/(pow((pow(x-earth_x,2)+pow(y-earth_y,2)),3/2))  -G * mm*(x-moon_x)/(pow(pow(x-moon_x,2)+pow(y-moon_y,2),3/2))
    y_acc= -G * me * (y-earth_y)/(pow((pow(x-earth_x,2)+pow(y-earth_y,2)),3/2))  -G * mm*(y-moon_y)/(pow((pow(x-moon_x,2)+pow(y-moon_y,2)),3/2))
    #print(x_acc,y_acc)
    #print(moon_x,t)
    #plt.scatter(earth_x,earth_y,color='red')
    #plt.scatter(moon_x,moon_y,color='black')
    if x_acc>0:
        print(t)
    return x_acc,y_acc 

def taylor(x,y,x_vel,y_vel,x_acc,y_acc,a,t):
    x1=x+a*x_vel+pow(a,2)/2*x_acc
    #print(x)
    #print(x_acc)
    #print(x1)
    y1=y+a*y_vel+pow(a,2)/2*y_acc
    x1_vel=x_vel+a*x_acc
    y1_vel=y_vel+a*y_acc
    x1_acc,y1_acc=acceleration(x1,y1,t+a)
    return x1,y1,x1_vel,y1_vel,x1_acc,y1_acc
    
def run_taylor(start,nsteps,a,x,v_init):
    t=0
    x1= x #r12+delta_r-re 
    y1=0
    x1_vel=0
    y1_vel=v_init
    x1_acc,y1_acc=acceleration(x,0,t)
    xs=[x]
    ys=[0]
    for i in range(nsteps):
        t=a*(i+1)
        x1,y1,x1_vel,y1_vel,x1_acc,y1_acc=taylor(x1,y1,x1_vel,y1_vel,x1_acc,y1_acc,a,t)
        xs.append(x1)
        ys.append(y1)
    return xs,ys

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

def run_runge(start,nsteps,a,x,v_init):
    t=0
    x1=x
    y1=0
    x1_vel=0
    y1_vel=v_init
    #print(y1_vel)
    x1s=[x1]
    y1s=[y1]
    x1_acc,y1_acc=acceleration(x,0,t)
    x1s_acc=[x1_acc]
    y1s_acc=[y1_acc]
    for i in range(nsteps):
        t=a*i
        x1,y1,x1_vel,y1_vel,x1_acc,y1_acc = runge_kutta(x1,y1,x1_vel,y1_vel,x1_acc,y1_acc,a,t)
        #print(y1_vel)
        #plt.scatter(x1,y1)
        x1s.append(x1)
        y1s.append(y1)
        x1s_acc.append(x1_acc)
        y1s_acc.append(y1_acc)
        #print(x1,y1,x1_vel,y1_vel,x1_acc,y1_acc)
    plt.plot(x1s,y1s)
    #print(a*nsteps)
 #print(x1s)
    #print(x1s_acc)
    return x1s,y1s


def run2():
    for i in range(10):
        nsteps=12000
        a=200
        pert=1.05951+0.000001*i
        delta_r=r12*pow(mm/(3*me),1/3)*pert
        x=r12+delta_r-re
        r=run_runge(0,nsteps,a,x,v_init)
    plt.scatter(-re*np.cos(2*np.pi*nsteps*a/T),-re*np.sin(2*np.pi*nsteps*a/T))
    plt.scatter(rm*np.cos(2*np.pi*nsteps*a/T),rm*np.sin(2*np.pi*nsteps*a/T))
    return(a*nsteps)

def orbit_2():
    nsteps=24000
    a=200
    pert=1.059515200719
    delta_r=r12*pow(mm/(3*me),1/3)*pert
    x=r12+delta_r-re
    r=run_runge(0,nsteps,a,x,v_init)
    plt.scatter(-re*np.cos(2*np.pi*nsteps*a/T),-re*np.sin(2*np.pi*nsteps*a/T))
    plt.scatter(rm*np.cos(2*np.pi*nsteps*a/T),rm*np.sin(2*np.pi*nsteps*a/T))


pert=  1.059515200719  #1.05945925 
delta_r=r12*pow(mm/(3*me),1/3)*pert
x=r12+delta_r-re
a=200
nsteps= int(T/a)
xs,ys=run_runge(0,nsteps,a,x,v_init)
plt.scatter(-re*np.cos(2*np.pi*nsteps*a/T),-re*np.sin(2*np.pi*nsteps*a/T),label="Earth",color='blue')
plt.scatter(rm*np.cos(2*np.pi*nsteps*a/T),rm*np.sin(2*np.pi*nsteps*a/T),label="Moon",color='black')
#plt.plot(delta_r,0,color='blue') 
plt.plot(xs,ys,color="red",linewidth=1.75,label='Rocket Orbit')
plt.xlabel('X Positions (m)')
plt.ylabel('Y Positions (m)')
plt.legend(loc='upper right')
plt.title('Rocket Trajectory Using RK4')
#plt.savefig('/Users/mac/Desktop/computing/final orbit.png',transparent=True)

# Set the background color to be transparent
#   # 'none' makes the background transparent

# Save the plot with a transparent background


print(delta_r)
#t=orbit_2()

#r=orbit_2()
#print('hello')
plt.show()