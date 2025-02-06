import math
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric, CartesianRepresentation, EarthLocation, AltAz
import astropy.units as u

# Asteroid Orbital Parameters
epoch_of_osculation = Time("2024-01-09")
mean_anomaly = math.radians(175.9094985980576)
argument_of_perihelion = math.radians(48.48336)
long_of_ascending_node = math.radians(247.23814)
inclination = math.radians(16.65913)
eccentricity = 0.1115598
semimajor_axis = 5.23060638 * u.au

# Observer Location (Durham)
durham_location = EarthLocation(lat='54:46:01', lon='-1:34:24', height=119.5 * u.meter)

# Dates for Observation
date_strings = [
    '2022-10-10T20:58:22.383'
    # Add more dates here...
]

# Initialize lists to store coordinates
asteroid_coordinates = []

for date_str in date_strings:
    observing_time = Time(date_str)

    # Calculate the time since epoch in seconds
    time_since_epoch = (observing_time - epoch_of_osculation).sec

    # Calculate the eccentric anomaly (E) using Newton-Raphson iteration
    E = mean_anomaly
    tolerance = 1e-9
    while True:
        E_next = E - (E - eccentricity * math.sin(E) - mean_anomaly) / (1 - eccentricity * math.cos(E))
        if abs(E_next - E) < tolerance:
            break
        E = E_next

    # Calculate the true anomaly (ν)
    nu = 2 * math.atan2(math.sqrt(1 + eccentricity) * math.sin(E / 2), math.sqrt(1 - eccentricity) * math.cos(E / 2))

    # Calculate the distance from the focus (r)
    r = (semimajor_axis * (1 - eccentricity * math.cos(E))).to(u.au)

    # Calculate the position in the orbital plane (x, y)
    x = (r * (math.cos(long_of_ascending_node) * math.cos(nu + argument_of_perihelion) - math.sin(long_of_ascending_node) * math.sin(nu + argument_of_perihelion))).to(u.au)
    y = (r * (math.sin(long_of_ascending_node) * math.cos(nu + argument_of_perihelion) + math.cos(long_of_ascending_node) * math.sin(nu + argument_of_perihelion))).to(u.au)

    # Calculate the z-coordinate (assuming the orbit is in the ecliptic plane)
    z = 0 * u.au  # Z-coordinate in astronomical units (au)

    # Create a CartesianRepresentation
    asteroid_position = CartesianRepresentation(x, y, z)

    # Calculate the position of the asteroid in the Solar System barycentric frame
    with solar_system_ephemeris.set('builtin'):
        asteroid_position_barycentric = get_body_barycentric('sun', observing_time) - asteroid_position

    # Append the coordinates to the list
    asteroid_coordinates.append((asteroid_position_barycentric.x.value, asteroid_position_barycentric.y.value, asteroid_position_barycentric.z.value))
    print(f"Date {date_str}:")
    print("X Coordinate:", asteroid_position_barycentric.x)
    print("Y Coordinate:", asteroid_position_barycentric.y)
    print("Z Coordinate:", asteroid_position_barycentric.z)
    print("-------------------")

# Print the coordinates for all dates
for i, (x, y, z) in enumerate(asteroid_coordinates, start=1):
    print(f"Date {i}:")
    print("X Coordinate:", x)
    print("Y Coordinate:", y)
    print("Z Coordinate:", z)
    print("-------------------")
