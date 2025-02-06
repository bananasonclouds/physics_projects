from ephem import *
import time, math, subprocess, sys, os
import numpy as np

'''date_obs = sys.argv[1]
year = date_obs[0:4]
month = date_obs[5:7]
day = date_obs[8:10]
obstime = date_obs[11:]'''

#print(year,month,day,obstime)

elements = open('test.txt','r')
minor = EllipticalBody()
minor._epoch = '2000/01/01.5'   # epoch date for J2000.0
ele = np.genfromtxt('test',dtype=None)
minor._M = float(ele[2].decode())
minor._om = float(ele[3].decode())
minor._Om = float(ele[4].decode())
minor._inc = float(ele[5].decode())
minor._e = float(ele[6].decode())
minor._a = float(ele[7].decode())
minor._H = float(ele[8].decode())
minor._G = float(ele[9].decode())
minor._epoch_M = ele[0].decode()#word[4]+'/'+word[5]+'/'+word[6]


durham = Observer()
durham.epoch = '2000'
durham.long = '-1:34:24'
durham.lat = '54:46:01'
durham.temp = 10
durham.elev = 119.5
durham.date = year+'/'+month+'/'+day+ ' '+obstime
#print(durham.date)

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
#time                dist    RA        Dec     xA[au]  yA[au]  zA[au]   dsol[au] RAsol[au]  DecSol[au] xSun[au] ySun[au] zSun[au]
print (sys.argv[1], "%7.5f" % d_minor, "%7.5f" % (180./math.pi*ecl_long_minor), "%7.5f" % (180./math.pi*ecl_lat_minor), "%7.5f" % x_minor, "%7.5f" %  y_minor, "%7.5f" % z_minor,  "%7.5f" % d_sun, "%7.5f" % (180./math.pi*ecl_long_sun), "%7.5f" % (180./math.pi*ecl_lat_sun), "%7.5f" % x_sun, "%7.5f" %  y_sun, "%7.5f" % z_sun)
