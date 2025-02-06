from ephem import *
import time, math, subprocess, sys, os
import numpy as np

# Define a list of date strings
date_strings = [
'2022-08-10 00:00:00',
'2022-10-09 00:00:00',
'2022-12-08 00:00:00',
'2023-02-06 00:00:00',
'2023-04-07 00:00:00',
'2023-06-06 00:00:00',
'2023-08-05 00:00:00',
'2023-10-04 00:00:00',
'2023-12-03 00:00:00',
'2024-02-01 00:00:00',
'2024-04-01 00:00:00',
'2024-05-31 00:00:00',
'2024-07-30 00:00:00',
'2024-09-28 00:00:00',
'2024-11-27 00:00:00',
'2025-01-26 00:00:00',
'2025-03-27 00:00:00',
'2025-05-26 00:00:00',
'2025-07-25 00:00:00',
'2025-09-23 00:00:00',
'2025-11-22 00:00:00',
'2026-01-21 00:00:00',
'2026-03-22 00:00:00',
'2026-05-21 00:00:00',
'2026-07-20 00:00:00',
'2026-09-18 00:00:00',
'2026-11-17 00:00:00',
'2027-01-16 00:00:00',
'2027-03-17 00:00:00',
'2027-05-16 00:00:00',
'2027-07-15 00:00:00',
'2027-09-13 00:00:00',
'2027-11-12 00:00:00',
'2028-01-11 00:00:00',
'2028-03-11 00:00:00',
'2028-05-10 00:00:00',
'2028-07-09 00:00:00',
'2028-09-07 00:00:00',
'2028-11-06 00:00:00',
'2029-01-05 00:00:00',
'2029-03-06 00:00:00',
'2029-05-05 00:00:00',
'2029-07-04 00:00:00',
'2029-09-02 00:00:00',
'2029-11-01 00:00:00',
'2029-12-31 00:00:00',
'2030-03-01 00:00:00',
'2030-04-30 00:00:00',
'2030-06-29 00:00:00',
'2030-08-28 00:00:00',
'2030-10-27 00:00:00',
'2030-12-26 00:00:00',
'2031-02-24 00:00:00',
'2031-04-25 00:00:00',
'2031-06-24 00:00:00',
'2031-08-23 00:00:00',
'2031-10-22 00:00:00',
'2031-12-21 00:00:00',
'2032-02-19 00:00:00',
'2032-04-19 00:00:00',
'2032-06-18 00:00:00',
'2032-08-17 00:00:00',
'2032-10-16 00:00:00',
'2032-12-15 00:00:00',
'2033-02-13 00:00:00',
'2033-04-14 00:00:00',
'2033-06-13 00:00:00',
'2033-08-12 00:00:00',
'2033-10-11 00:00:00',
'2033-12-10 00:00:00',
'2034-02-08 00:00:00',
'2034-04-09 00:00:00',
'2034-06-08 00:00:00',
'2034-08-07 00:00:00',
'2034-10-06 00:00:00',
'2034-12-05 00:00:00',
'2035-02-03 00:00:00',
'2035-04-04 00:00:00',
'2035-06-03 00:00:00',
'2035-08-02 00:00:00',
'2035-10-01 00:00:00',
'2035-11-30 00:00:00',
'2036-01-29 00:00:00',
'2036-03-29 00:00:00',
'2036-05-28 00:00:00',
'2036-07-27 00:00:00',
'2036-09-25 00:00:00',
'2036-11-24 00:00:00',
'2037-01-23 00:00:00',
'2037-03-24 00:00:00',
'2037-05-23 00:00:00',
'2037-07-22 00:00:00',
'2037-09-20 00:00:00',
'2037-11-19 00:00:00',
'2038-01-18 00:00:00',
'2038-03-19 00:00:00',
'2038-05-18 00:00:00',
'2038-07-17 00:00:00',
'2038-09-15 00:00:00',
'2038-11-14 00:00:00'

  # Add more dates here...
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

    elements = open('/Users/mac/Desktop/ORBITAL P/test_kallisto','r')
    minor = EllipticalBody()
    minor._epoch = '2000/01/01.5'   # epoch date for J2000.0
    ele = np.genfromtxt('/Users/mac/Desktop/ORBITAL P/test_kallisto',dtype=None)
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
    #print("Date:", date)
    #print("Distance to Minor Body:", d_minor)
    #print("Ecliptic Longitude of Minor Body:", ecl_long_minor)
    #print("Ecliptic Latitude of Minor Body:", ecl_lat_minor)
    print("X_Minor Body:", x_minor)
    print("Y_Minor Body:", y_minor)
    print("Z_Minor Body:", z_minor)
    #print("Distance to Sun:", d_sun)
    #print("Ecliptic Longitude of Sun:", ecl_long_sun)
    #print("Ecliptic Latitude of Sun:", ecl_lat_sun)
    print("X_Sun:", x_sun)
    print("Y_Sun:", y_sun)
    print("Z_Sun:", z_sun)
    print("-------------------")

