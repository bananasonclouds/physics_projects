def getrot(hdr):

    cd11 = hdr['CD1_1']
    cd21 = hdr['CD2_1']

    rho_a = math.atan2(cd21, cd11)
    angle = 180-np.rad2deg(rho_a)
    return angle

import astropy
import numpy as np
import math
from astropy.io import fits

hdu = fits.open('/mnt/archive/draco2/2023/23_11_28/ad0008.fits')
#hdu = fits.open('/mnt/archive/draco2/2023/23_01_07/ad0035.fits')
img = hdu[0].data
hdr = hdu[0].header

print(getrot(hdr))

rotation = getrot(hdr)

def rotate(origin, point, angle):

    rad_angle = np.radians(angle)

    ox, oy = origin
    px, py = point

    qx = ox + np.cos(rad_angle) * (px-ox) - np.sin(rad_angle) * (py - oy)
    qy = oy + np.sin(rad_angle) * (px - ox) + np.cos(rad_angle) *(py-oy)

    qx = qx - ox
    qy = qy - oy

    return qx, qy

uranus = 758.65, 576.35
#miranda = 551.1, 410.0
oberon = 748.61, 531.87
titania = 739.9, 550.34
#ariel = 635.45293, 296.29745
umbriel = 771.1, 594.3

#print(rotate(uranus, miranda, rotation))
print(rotate(uranus, oberon, rotation))
print(rotate(uranus, titania, rotation))
#print(rotate(uranus, ariel , rotation))
print(rotate(uranus, umbriel, rotation))

neptune = 562.25, 429.25
triton = 570.64, 427.5

print(rotate(neptune, triton, rotation))



