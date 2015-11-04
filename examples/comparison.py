# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import subprocess
import struct

import gridmap
try:
    import pyproj
    has_pyproj = True
except ImportError:
    print('pyproj is not installed')
    has_pyproj = False
# from mpl_toolkits.basemap import Basemap

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 1000        # grid resolution (at lat_ts)  [m]
ylon  = 58.0          # angle of y-axis        [deg]

x, y = 200.0, 100.0

print("\n --- Spherical earth ---\n")

# Define objects

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)

if has_pyproj:
    pmap = pyproj.Proj(gmap.proj4string)

# bmap = Basemap(projection='stere', rsphere=6371000,
#                llcrnrlon = ylon, llcrnrlat = 90.0,
#                urcrnrlon = ylon+1, urcrnrlat = 89.0,
#                lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)

def proj(projstring, lon, lat):
    """Use proj command for forward projection"""
    command = 'proj -o ' + projstring
    p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE,
                     stdout=subprocess.PIPE)
    p.stdin.write("%s %s\n" % (str(lon), str(lat)))
    out, err = p.communicate()
    return struct.unpack('2d', out)

def gmt_forward(gmap, lon, lat):
    """Use GMT mapproject for forward projection"""

    # m = 193.713819606 / 208.754891658  # Hva fanden?
    m = 1
    # m = 1. / (1 +np.cos(np.pi / 3))  # 60 degrees
    # m = 2./  (1 +np.cos(np.pi / 3))  # 60 degrees
    # m = 309.672008317 / 283.177139485
    projection = '-Js%s/90.0/%s/1:%s -Dm' % \
                 (str(ylon), str(gmap.lat_ts), str(m*dx))
    # projection = '-Js%s/90.0/%s/1:%s -F' % \
    #         (str(ylon), str(gmap.lat_ts), str(2*dx))

    #ellipsoid = ""
    #ellipsoid = "--ELLIPSOID=Sphere"
    #ellipsoid = "--ELLIPSOID=WGS-84"
    ## Ser ut som ellipsoid er gitt på korrekt vis, men får feil
    ## Setting virker for WGS84
    #ellipsoid = "--ELLIPSOID=%s" % str(gmap.ellipsoid.a)
    #ellipsoid = "--ELLIPSOID=%s" % "6378137,298.257223563" # WGS84 OK
    # Den under virker og gir WGS84
    #ellipsoid = "--ELLIPSOID=%s%s" % ("6378137,f=", str(1/298.257223563))
    #ellipsoid = "--ELLIPSOID=%s" % "6371000,f=0"
    ellipsoid = "--ELLIPSOID=%s" % "6371000,b=6371000"

    m = 1 / (1 + np.sin(np.pi / 3))
    # extent = '-R0/1/60/61 --MAP_SCALE_FACTOR=%s' % str(m)   # Actual values are not used
    # extent = '-R0/1/60/61'   # Actual values are not used
    extent = '-Rg'
    # offset = '-C%s/%s' % (str(dx*xp), str(dx*yp))
    offset=" "
    command0 = 'GMT mapproject -bo'
    gmtstring = " ".join((ellipsoid, projection, offset, extent))
    command = " ".join((command0, ellipsoid, projection, offset, extent))
    # Set up the  process
    # if verbose: print(command)
    print(command)
    p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE,
                     stdout=subprocess.PIPE)

    # Send lon, lat to mapproject
    # p.stdin.write("{:f} {:f}\n".format(lon, lat))
    p.stdin.write("%s %s\n" % (str(lon), str(lat)))


    # Get the output
    out, err = p.communicate()
    # print("err, out = ", err, out)
    x, y = struct.unpack('dd', out)
    #x = x / (m*dx) + xp
    #y = y / (m*dx) + yp
    return x, y


# --- Forward projection

lon, lat = 2, 66   # Station M
lon, lat = 0, 90   # Station M

print("gmap.ll2grid   : ", *gmap.ll2grid(lon, lat))

if has_pyproj:
    px, py = pmap(lon, lat)
    print("pyproj         : ", px/dx, py/dx)

# basemap
# bx, by = bmap(lon, lat)
# print("basemap        : ", bx/dx + xp, by/dx + yp)

# proj4
px, py = proj(gmap.proj4string, lon, lat)
print("proj4          : ", px/dx, py/dx)
#print("difference [m] : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * gmap.dx)

# GMT
print("GMT mapproject : ", *gmt_forward(gmap, lon, lat))


import sys
sys.exit()


print("\n --- WGS84 ---\n")

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, ellipsoid='WGS84')

pmap = pyproj.Proj(gmap.proj4string)

x0, y0 = gmap.ll2grid(lon, lat)
x1, y1 = pmap(lon, lat)
x1, y1 = x1/gmap.dx, y1/gmap.dx

print("pyproj           : ", x1, y1)
print("gmap.ll2grid   : ", x0, y0)
print("difference [m] : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * gmap.dx)


lon0, lat0 = gmap.grid2ll(x, y)
lon1, lat1 = pmap(x*gmap.dx, y*gmap.dx, inverse=True)
print()
print("pyproj         : ", lon1, lat1)
print("gmap.grid2ll   : ", lon0, lat0)

lon1, lat1 = pmap(x*gmap.dx, y*gmap.dx, inverse=True)
print()
print("invproj        : ", lon1, lat1)
print("gmap.grid2ll   : ", lon0, lat0)
