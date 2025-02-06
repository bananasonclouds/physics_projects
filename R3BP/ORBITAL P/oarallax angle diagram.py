import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Create figure and axis
fig, ax = plt.subplots()

# Define coordinates and sizes
star_pos = (2, 4)
baseline_half_length = 1

# Draw Star
star = plt.Circle(star_pos, 0.1, color='yellow', label='Distant Star')
ax.add_patch(star)

# Draw baseline
baseline_start = (-baseline_half_length, 0)
baseline_end = (baseline_half_length, 0)
ax.plot([baseline_start[0], baseline_end[0]], [baseline_start[1], baseline_end[1]], color='green', label='Baseline')

# Draw sight lines (representing the hypotenuse)
ax.plot([baseline_start[0], star_pos[0]], [baseline_start[1], star_pos[1]], linestyle='--', color='red')
ax.plot([baseline_end[0], star_pos[0]], [baseline_end[1], star_pos[1]], linestyle='--', color='red')

# Label points and distances
ax.text(star_pos[0], star_pos[1] + 0.2, 'Star', horizontalalignment='center')
ax.text(baseline_start[0], baseline_start[1] - 0.2, 'Durham', horizontalalignment='right')
ax.text(baseline_end[0], baseline_end[1] - 0.2, 'La Palma', horizontalalignment='left')
ax.text(star_pos[-2]-2, star_pos[0] - 1.5, r'$\alpha$', color='orange')

# Add angle arc for parallax
parallax_arc = patches.Arc((-1, 0), 2 * baseline_half_length, 2 * baseline_half_length,
                           angle=0, theta1=0, theta2=45, edgecolor='orange', label='Parallax Angle')
ax.add_patch(parallax_arc)

# Set limits and aspect
ax.set_xlim(-2, 4)
ax.set_ylim(-1, 5)
ax.set_aspect('equal')

# Show the plot
plt.title("Parallax Triangulation Method")
plt.show()

import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Create figure and axis
fig, ax = plt.subplots()

# Define coordinates and sizes
star_pos = (2, 4)
baseline_half_length = 1

# Draw Star
star = plt.Circle(star_pos, 0.1, color='yellow', label='Distant Star')
ax.add_patch(star)

# Draw baseline
baseline_start = (-baseline_half_length, 0)
baseline_end = (baseline_half_length, 0)
ax.plot([baseline_start[0], baseline_end[0]], [baseline_start[1], baseline_end[1]], color='green', label='Baseline')

# Draw sight lines (representing the hypotenuse)
ax.plot([baseline_start[0], star_pos[0]], [baseline_start[1], star_pos[1]], linestyle='--', color='red')
ax.plot([baseline_end[0], star_pos[0]], [baseline_end[1], star_pos[1]], linestyle='--', color='red')

# Label points and distances
ax.text(star_pos[0], star_pos[1] + 0.2, 'Star', horizontalalignment='center')
ax.text(baseline_start[0], baseline_start[1] - 0.2, 'Durham', horizontalalignment='right')
ax.text(baseline_end[0], baseline_end[1] - 0.2, 'La Palma', horizontalalignment='left')

# Add angle arc for parallax
parallax_arc = patches.Arc((0, 0), 2 * baseline_half_length, 2 * baseline_half_length,
                           angle=0, theta1=0, theta2=45, edgecolor='orange', label='Parallax Angle')
ax.add_patch(parallax_arc)

# Add line for right angle
right_angle_line = patches.ConnectionPatch(xyA=baseline_start, xyB=star_pos, coordsA="data", coordsB="data",
                                           arrowstyle='-', linestyle='dotted', color='black')
ax.add_patch(right_angle_line)

# Set limits and aspect
ax.set_xlim(-2, 4)
ax.set_ylim(-1, 5)
ax.set_aspect('equal')

# Show the plot
plt.title("Parallax Triangulation Method")
plt.show()
