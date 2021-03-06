# -*- coding: utf-8 -*-

from __future__ import print_function

import sys
import subprocess
import struct
from math import pi
#sys.path = ['../gridmap'] + sys.path # import from developing version
from gridmap import *

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_c)   [m]
ylon  = 58.0          # angle of y-axis        [deg]

lon, lat = 2, 66   # Station M
x, y = 200.0, 100.0

verbose = True

# -----------------------------------------------------

def proj(proj4string, lon, lat):
    """Use proj4 for forward projection"""

    # Use -o for binary outpur => no errors due to format
    command = 'proj -o ' + proj4string

    # Set up the proj process
    if verbose: print(command)
    p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

    # Send lon, lat to proj
    p.stdin.write("%s %s\n" % (str(lon), str(lat)))

    # Get the output 
    out, err = p.communicate()
    x, y = struct.unpack('2d', out)
    return x, y


def invproj(proj4string, x, y):
    """Use proj4 for inverse projection"""

    command = 'invproj -o ' + proj4string

    if verbose: print(command)
    p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

    p.stdin.write("%s %s\n" % (str(x), str(y)))

    out, err = p.communicate()
    lon, lat = struct.unpack('2d', out)
    # Convert from radians to degrees
    lon = lon * 180.0 / pi
    lat = lat * 180.0 / pi

    return lon, lat

# ----------------------------------------------------

print("\n --- sphere ---\n")

gmap = PolarStereographic(xp, yp, dx, ylon)

x0, y0 = gmap.ll2grid(lon, lat)
x2, y2 = proj(gmap.proj4string2(), lon, lat)
x2, y2 = x2/gmap.dx, y2/gmap.dx


print("gmap.ll2grid : ", x0, y0)
print("proj4string1  : ", x1, y1)
print("proj4string2  : ", x2, y2)

print("distance in meter : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * dx)
print("distance in meter : ", ((x2-x0)**2 + (y2-y0)**2)**0.5 * dx)


print("\n --- WGS84 ---\n")

gmap = PolarStereographic(xp, yp, dx, ylon, ellipsoid=WGS84)

x0, y0 = gmap.ll2grid(lon, lat)
x1, y1 = proj(gmap.proj4string1(), lon, lat)
x2, y2 = proj(gmap.proj4string2(), lon, lat)
x2, y2 = x2/gmap.dx, y2/gmap.dx


print("gmap.ll2grid : ", x0, y0)
print("proj4string2  : ", x2, y2)

print("distance in meter : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * dx)
print("distance in meter : ", ((x2-x0)**2 + (y2-y0)**2)**0.5 * dx)

print("\n --- inverse sphere ---\n")

gmap = PolarStereographic(xp, yp, dx, ylon)

lon0, lat0 = gmap.grid2ll(x, y)
lon1, lat1 = invproj(gmap.proj4string1(), x, y)
lon2, lat2 = invproj(gmap.proj4string2(), x*dx, y*dx)

print("gmap.ll2grid : ", lon0, lat0)
print("proj4string1  : ", lon1, lat1)
print("proj4string2  : ", lon2, lat2)


print("\n --- inverse WGS84 ---\n")

gmap = PolarStereographic(xp, yp, dx, ylon, ellipsoid=WGS84)

lon0, lat0 = gmap.grid2ll(x, y)
lon1, lat1 = invproj(gmap.proj4string1(), x, y)
lon2, lat2 = invproj(gmap.proj4string2(), x*dx, y*dx)

print("gmap.ll2grid : ", lon0, lat0)
print("proj4string1  : ", lon1, lat1)
print("proj4string2  : ", lon2, lat2)






