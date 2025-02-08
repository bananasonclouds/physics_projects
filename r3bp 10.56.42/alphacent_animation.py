import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

# Define universal gravitation constant
G = 6.67408e-11 # N-m2/kg2

# Mass of the sun, and other masses
m_nd = 1.989e30  # Sun
m1 = 1.1         # Rigil Kentaurus
m2 = 0.907       # Toliman
m3 = 0.256       # Proxima Centauri

# Distance and orbital parameters
r_nd = 5.326e12 #distance between stars 
t_nd = 79*365*24*3600 #orbital period
v_nd = 23000 #orbital velocity

# Net constants
K1 = G*t_nd*m1*m_nd/((r_nd**2)*v_nd)
K2 = v_nd*t_nd/r_nd

# Define initial position vectors
def roll():
    return 2 * np.random.random_sample(3) - 1

r1, r2, r3 = roll(), roll(), roll()

# Define initial velocities
v1 = np.array([0.01, 0.01, 0])
v2 = np.array([-0.05, 0, -0.01])
v3 = np.array([0, -0.01, 0])

# Centre of Mass calculations
r_com = (m1*r1 + m2*r2 + m3*r3) / (m1 + m2 + m3)
v_com = (m1*v1 + m2*v2 + m3*v3) / (m1 + m2 + m3)

def ThreeBodyEquations(w, t, G, m1, m2, m3):
    DIMS = 3

    # We can insert a 0 value at "V0" for readability.
    # This is done to access v values with a
    # 1-based index only and is not a good practice.

    v = [0]

    for i in range(DIMS*2):
      v.append(w[i*DIMS:(i+1)*DIMS])

    norm = np.linalg.norm
    r12 = norm(v[2] - v[1])
    r13 = norm(v[3] - v[1])
    r23 = norm(v[3] - v[2])

    def dvdt(K1, m1, m2, v1, v2, v3, r1, r2):
      A = m1*(v2-v1)/r1**3 
      B = m2*(v3-v1)/r2**3
      C = K1
      return C*(A+B)

    dv1dt = dvdt(K1, m2, m3, v[1], v[2], v[3], r12, r13)
    dv2dt = dvdt(K1, m1, m3, v[2], v[1], v[3], r12, r23)
    dv3dt = dvdt(K1, m1, m2, v[3], v[1], v[2], r13, r23)

    dr1dt=K2*v[4]
    dr2dt=K2*v[5]
    dr3dt=K2*v[6]

    concat = np.concatenate

    r12_d = concat((dr1dt, dr2dt))
    r_d = concat((r12_d, dr3dt))
    v12_d = concat((dv1dt, dv2dt))
    v_d = concat((v12_d, dv3dt))
    d = concat((r_d, v_d))

    return d

# Package initial parameters and solve ODE
init_params = np.array([r1, r2, r3, v1, v2, v3]).flatten()
time_span = np.linspace(0, 250, 20000) 
three_body_sol = scipy.integrate.odeint(ThreeBodyEquations, init_params, time_span, args=(G, m1, m2, m3))

def extract_coords(r):
    l = []
    for i in range(3):
     l.append(r[:, i])
    return np.array(l)

coords_1 = extract_coords(three_body_sol[:, :3])
coords_2 = extract_coords(three_body_sol[:, 3:6])
coords_3 = extract_coords(three_body_sol[:, 6:9])

def update_lines(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines

# Plotting and animation setup
with plt.figure() as fig:
    ax = p3.Axes3D(fig)

# Data with position coordinates:
data = [coords_1, coords_2, coords_3]

# Creating line objects:
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

# Setting the axes properties
ax_limits = [-1, 1]
ax.set_xlim3d(ax_limits)
ax.set_xlabel('X')

ax.set_ylim3d(ax_limits)
ax.set_ylabel('Y')

ax.set_zlim3d(ax_limits)
ax.set_zlabel('Z')

ax.set_title('3D 3-Body Problem')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig,
                                   update_lines,
                                   fargs=(data, lines),
                                   frames = 500,
                                   interval = 10,
                                   blit=False,
)

# Display the animation (only works in IPython environments)
try:
    from IPython.display import HTML
    HTML(line_ani.to_html5_video())
except ImportError:
    print("IPython display is not available.")
