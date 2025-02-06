from ephem import *
import time, math, subprocess, sys, os
import numpy as np

# Define a list of date strings
date_strings = [
'2030-09-27 00:00:00'  # Add more dates here...
]

coordinates_list=[]

for date_str in date_strings:
    # Create an Observer for Durham with the specified date
    durham = Observer()
    durham.epoch = '2000'
    durham.long = '-1:34:24'
    durham.lat = '54:46:01'
    durham.temp = 10
    durham.elev = 119.5
    durham.date = date_str

    elements = open('/Users/mac/Desktop/ORBITAL P/test_helio','r')
    minor = EllipticalBody()
    minor._epoch = '2000/01/01.5'   # epoch date for J2000.0
    ele = np.genfromtxt('/Users/mac/Desktop/ORBITAL P/test_helio',dtype=None)
    minor._M = float(ele[2].decode())
    minor._om = float(ele[3].decode())
    minor._Om = float(ele[4].decode())
    minor._inc = float(ele[5].decode())
    minor._e = float(ele[6].decode())
    minor._a = float(ele[7].decode())
    minor._H = float(ele[8].decode())
    minor._G = float(ele[9].decode())
    minor._epoch_M = ele[0].decode()#word[4]+'/'+word[5]+'/'+word[6]

    # Perform the calculations for the minor body and the Sun as before
    const = 180./math.pi
    minor.compute(durham)
    pos = Equatorial(minor.a_ra,minor.a_dec, epoch='2000')
    #print 180./math.pi*float(minor.a_ra),  180./math.pi*float(minor.a_dec)
    ecl = Ecliptic(pos)
    ecl_long_minor = float(ecl.lon)
    ecl_lat_minor  = float(ecl.lat)

    d_minor = minor.earth_distance
    x_minor = d_minor*math.cos(ecl_lat_minor)*math.cos(ecl_long_minor)
    y_minor = d_minor*math.cos(ecl_lat_minor)*math.sin(ecl_long_minor)
    z_minor = d_minor*math.sin(ecl_lat_minor)

    xxx = Sun()
    xxx.compute(durham)
    pos = Equatorial(xxx.a_ra,xxx.a_dec, epoch='2000')
    ecl = Ecliptic(pos)
    ecl_long_sun = float(ecl.lon)
    ecl_lat_sun  = float(ecl.lat)
    #print ecl_long*180./math.pi , ecl_lat*180./math.pi 
    d_sun = xxx.earth_distance
    x_sun = d_sun*math.cos(ecl_lat_sun)*math.cos(ecl_long_sun)
    y_sun = d_sun*math.cos(ecl_lat_sun)*math.sin(ecl_long_sun)
    z_sun = d_sun*math.sin(ecl_lat_sun)

    
    coordinates_list.append((durham.date, d_minor, ecl_long_minor, ecl_lat_minor, x_minor, y_minor, z_minor, d_sun, ecl_long_sun, ecl_lat_sun, x_sun, y_sun, z_sun))

# Print the coordinates for all dates
for i, (date, d_minor, ecl_long_minor, ecl_lat_minor, x_minor, y_minor, z_minor, d_sun, ecl_long_sun, ecl_lat_sun, x_sun, y_sun, z_sun) in enumerate(coordinates_list, start=1):
    print(f"Date {i}:")
    print("Date:", date)
    print("Distance to Minor Body:", d_minor)
    print("Ecliptic Longitude of Minor Body:", ecl_long_minor)
    print("Ecliptic Latitude of Minor Body:", ecl_lat_minor)
    print("X_Minor Body:", x_minor)
    print("Y_Minor Body:", y_minor)
    print("Z_Minor Body:", z_minor)
    print("Distance to Sun:", d_sun)
    print("Ecliptic Longitude of Sun:", ecl_long_sun)
    print("Ecliptic Latitude of Sun:", ecl_lat_sun)
    print("X_Sun:", x_sun)
    print("Y_Sun:", y_sun)
    print("Z_Sun:", z_sun)
    print("-------------------")

