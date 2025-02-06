import numpy as np
import ephem
import math

def vary_orbital_parameters(parameters, errors, num_simulations):
    all_simulations = []
    for _ in range(num_simulations):
        varied_params = {}
        for param, value in parameters.items():
            error = errors.get(param, 0)
            varied_value = value + np.random.normal(0, scale=1) * error
            varied_params[param] = varied_value
        all_simulations.append(varied_params)
    return all_simulations


date_strings = [
'2022-08-10 00:00:00',
]

#aneas
orbital_parameters = {'semi_major_axis': 5.23060638 , 'eccentricity': 0.1115598, 'inclination': 16.65913, 'long': 247.23814, 'arg': 48.48336, 'MA': 49.56596635, 'Mag': 8.68, 'G': 0.15}
errors = {'semi_major_axis': 0.000800488, 'eccentricity': 0.006098104, 'inclination': 0.001550846, 'long': 0.031273521, 'arg': 2.127861057, 'MA': 1.338881049, 'Mag': 0, 'G': 0}
	
						
num_simulations = 300
varied_parameters = vary_orbital_parameters(orbital_parameters, errors, num_simulations)

# Iterate over each simulation and calculate coordinates
for i, simulation in enumerate(varied_parameters):
    print(f"Simulation {i+1}:")
    for date_str in date_strings:
        durham = ephem.Observer()
        durham.epoch = '2000'
        durham.long = '-1:34:24'
        durham.lat = '54:46:01'
        durham.temp = 10
        durham.elev = 119.5
        durham.date = date_str
        
        # Set parameters from the simulation
        minor = ephem.EllipticalBody()
        minor._epoch = '2000/01/01.5'
        minor._M = simulation['MA']
        minor._Om = simulation['long']
        minor._inc = simulation['inclination']
        minor._e = simulation['eccentricity']
        minor._a = simulation['semi_major_axis']
        minor._H = simulation['Mag']
        minor._G = simulation['G']
        minor._epoch_M = '2000/01/01.5'

        # Perform calculations
        minor.compute(durham)
        pos = ephem.Equatorial(minor.a_ra, minor.a_dec, epoch='2000')
        ecl = ephem.Ecliptic(pos)
        ecl_long_minor = float(ecl.lon)
        ecl_lat_minor = float(ecl.lat)

        d_minor = minor.earth_distance
        x_minor = d_minor * math.cos(ecl_lat_minor) * math.cos(ecl_long_minor)
        y_minor = d_minor * math.cos(ecl_lat_minor) * math.sin(ecl_long_minor)
        z_minor = d_minor * math.sin(ecl_lat_minor)

        print(f"Date: {date_str}, X: {x_minor}, Y: {y_minor}, Z: {z_minor}")
    print("-------------------")
