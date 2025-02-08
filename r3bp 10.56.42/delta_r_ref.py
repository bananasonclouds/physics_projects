import numpy as np
import matplotlib.pyplot as plt

def calculate_rocket_moon_separation():
    me = 5.9742 * 10 ** 24  # Mass of the Earth in kg
    mm = 7.35 * 10 ** 22  # Mass of the Moon in kg
    G = 6.6726 * 10 ** -11  # Gravitational constant in m^3 kg^-1 s^-2
    r12 = 384400.0 * 10 ** 3  # Distance between the CoM of the Earth and Moon in meters
    M = me + mm

    re = r12 * mm / (mm + me)
    rm = r12 * me / (mm + me)

    m_tot = me + mm
    delta_r = r12 * pow(mm / (3 * me), 1 / 3)
    T = np.sqrt(4 * pow(np.pi, 2) * pow(r12, 3) / (G * m_tot))

    def acceleration(x, y, t):
        earth_x = -re * np.cos(2 * np.pi * t / T)
        earth_y = -re * np.sin(2 * np.pi * t / T)
        moon_x = rm * np.cos(2 * np.pi * t / T)
        moon_y = rm * np.sin(2 * np.pi * t / T)
        x_acc = -G * me * (x - earth_x) / (pow((pow(x - earth_x, 2) + pow(y - earth_y, 2)), 3 / 2)) - G * mm * (
                    x - moon_x) / (pow(pow(x - moon_x, 2) + pow(y - moon_y, 2), 3 / 2))
        y_acc = -G * me * (y - earth_y) / (pow((pow(x - earth_x, 2) + pow(y - earth_y, 2)), 3 / 2)) - G * mm * (
                    y - moon_y) / (pow((pow(x - moon_x, 2) + pow(y - moon_y, 2)), 3 / 2))
        return x_acc, y_acc

    def taylor_expansion(x, y, x_vel, y_vel, a, t):
        x1 = x + a * x_vel + 0.5 * a**2 * acceleration(x, y, t)[0]
        y1 = y + a * y_vel + 0.5 * a**2 * acceleration(x, y, t)[1]
        x1_vel = x_vel + a * acceleration(x, y, t)[0]
        y1_vel = y_vel + a * acceleration(x, y, t)[1]
        return x1, y1, x1_vel, y1_vel

    def runge_kutta(x, y, x_vel, y_vel, a, t):
        x1 = x + a * x_vel
        y1 = y + a * y_vel
        x1_vel, y1_vel = acceleration(x1, y1, t + a / 2)
        return x1, y1, x1_vel, y1_vel

    t = 0
    x = re + delta_r  # Initial rocket position
    y = 0
    x_vel = 0
    y_vel = 2 * np.pi * (r12 + delta_r - re) / T  # Initial rocket velocity

    separation_data = []  # To store separation data
    time_data = []  # To store time data

    delta_t_values = [1, 10, 100, 1000]  # Different time step values to test

    for delta_t in delta_t_values:
        t = 0
        x = re + delta_r
        y = 0
        x_vel = 0
        y_vel = 2 * np.pi * (r12 + delta_r - re) / T
        separation_data = []  # Reset separation data for each delta_t
        time_data = []  # Reset time data for each delta_t

        while t < 27.3 * 24 * 3600:  # Stop after one orbit (27.3 days)
            x, y, x_vel, y_vel = taylor_expansion(x, y, x_vel, y_vel, delta_t, t)
            separation = np.sqrt((x - rm * np.cos(2 * np.pi * t / T)) ** 2 + (y - rm * np.sin(2 * np.pi * t / T)) ** 2)
            separation_data.append(separation)
            time_data.append(t)
            t += delta_t

        # Calculate L2_ref for the given delta_t
        L2_ref = separation_data[-1]  # Use the final separation as L2_ref

        # Plot (r_sat - L2_ref) / r_sat as a function of time for the current delta_t
        plt.plot(time_data, [(r - L2_ref) / r for r in separation_data], label=f'Delta_t = {delta_t} s')

    plt.xlabel('Time (s)')
    plt.ylabel('(r_sat - L2_ref) / r_sat')
    plt.title('Convergence for Different Time Steps')

    plt.legend()
    plt.show()
    plt.savefig('/Users/mac/Desktop/computing/final r_ref.png',transparent=True)


# Call the function to calculate and plot the convergence
calculate_rocket_moon_separation()
