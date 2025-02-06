import math
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric, CartesianRepresentation, EarthLocation, AltAz
import astropy.units as u

#Aneas Orbital Parameters
epoch_of_osculation = Time("2024-01-09")
mean_anomaly = math.radians(175.9094985980576)
argument_of_perihelion = math.radians(48.48336)
long_of_ascending_node = math.radians(247.23814)
inclination = math.radians(16.65913)
eccentricity = 0.1115598
semimajor_axis = 5.23060638 * u.au 

durham_location = EarthLocation(
    lat='54:46:01', lon='-1:34:24', height=119.5 * u.meter)

#REMINDER:delete last comma

date_strings = [
'2022-08-10T00:00:00',
'2022-09-09T00:00:00',
'2022-10-09T00:00:00',
'2022-11-08T00:00:00',
'2022-12-08T00:00:00',
'2023-01-07T00:00:00',
'2023-02-06T00:00:00',
'2023-03-08T00:00:00',
'2023-04-07T00:00:00',
'2023-05-07T00:00:00',
'2023-06-06T00:00:00',
'2023-07-06T00:00:00',
'2023-08-05T00:00:00',
'2023-09-04T00:00:00',
'2023-10-04T00:00:00',
'2023-11-03T00:00:00',
'2023-12-03T00:00:00',
'2024-01-02T00:00:00',
'2024-02-01T00:00:00',
'2024-03-02T00:00:00',
'2024-04-01T00:00:00',
'2024-05-01T00:00:00',
'2024-05-31T00:00:00',
'2024-06-30T00:00:00',
'2024-07-30T00:00:00',
'2024-08-29T00:00:00',
'2024-09-28T00:00:00',
'2024-10-28T00:00:00',
'2024-11-27T00:00:00',
'2024-12-27T00:00:00',
'2025-01-26T00:00:00',
'2025-02-25T00:00:00',
'2025-03-27T00:00:00',
'2025-04-26T00:00:00',
'2025-05-26T00:00:00',
'2025-06-25T00:00:00',
'2025-07-25T00:00:00',
'2025-08-24T00:00:00',
'2025-09-23T00:00:00',
'2025-10-23T00:00:00',
'2025-11-22T00:00:00',
'2025-12-22T00:00:00',
'2026-01-21T00:00:00',
'2026-02-20T00:00:00',
'2026-03-22T00:00:00',
'2026-04-21T00:00:00',
'2026-05-21T00:00:00',
'2026-06-20T00:00:00',
'2026-07-20T00:00:00',
'2026-08-19T00:00:00',
'2026-09-18T00:00:00',
'2026-10-18T00:00:00',
'2026-11-17T00:00:00',
'2026-12-17T00:00:00',
'2027-01-16T00:00:00',
'2027-02-15T00:00:00',
'2027-03-17T00:00:00',
'2027-04-16T00:00:00',
'2027-05-16T00:00:00',
'2027-06-15T00:00:00',
'2027-07-15T00:00:00',
'2027-08-14T00:00:00',
'2027-09-13T00:00:00',
'2027-10-13T00:00:00',
'2027-11-12T00:00:00',
'2027-12-12T00:00:00',
'2028-01-11T00:00:00',
'2028-02-10T00:00:00',
'2028-03-11T00:00:00',
'2028-04-10T00:00:00',
'2028-05-10T00:00:00',
'2028-06-09T00:00:00',
'2028-07-09T00:00:00',
'2028-08-08T00:00:00',
'2028-09-07T00:00:00',
'2028-10-07T00:00:00',
'2028-11-06T00:00:00',
'2028-12-06T00:00:00',
'2029-01-05T00:00:00',
'2029-02-04T00:00:00',
'2029-03-06T00:00:00',
'2029-04-05T00:00:00',
'2029-05-05T00:00:00',
'2029-06-04T00:00:00',
'2029-07-04T00:00:00',
'2029-08-03T00:00:00',
'2029-09-02T00:00:00',
'2029-10-02T00:00:00',
'2029-11-01T00:00:00',
'2029-12-01T00:00:00',
'2029-12-31T00:00:00',
'2030-01-30T00:00:00',
'2030-03-01T00:00:00',
'2030-03-31T00:00:00',
'2030-04-30T00:00:00',
'2030-05-30T00:00:00',
'2030-06-29T00:00:00',
'2030-07-29T00:00:00',
'2030-08-28T00:00:00',
'2030-09-27T00:00:00'
''   # Add more dates here...
]

coordinates_list=[]

for date_str in date_strings:
    observing_time = Time(date_str)

    durham = AltAz(location=durham_location, obstime=observing_time)

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
    position = CartesianRepresentation(x, y, z)

    # Calculate the position of the object in the Solar System barycentric frame
    with solar_system_ephemeris.set('builtin'):
        object_position = get_body_barycentric('earth', observing_time)  # Change 'earth' to the desired target body if needed

    # Add the position in the orbital plane to the Solar System barycentric position
    object_position += position

    # Append the coordinates to the list
    coordinates_list.append((object_position.x, object_position.y, object_position.z))
    print(f"Date {date_str}:")
    print("X Coordinate:", x)
    print("Y Coordinate:", y)
    print("Z Coordinate:", z)
    print("-------------------")

# Print the coordinates for all dates
for i, (x, y, z) in enumerate(coordinates_list, start=1):
    print(f"Date {i}:")
    print("X Coordinate:", x)
    print("Y Coordinate:", y)
    print("Z Coordinate:", z)
    print("-------------------")

