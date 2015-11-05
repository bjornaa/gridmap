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
dx    = 10000        # grid resolution (at lat_ts)  [m]
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

# TODO: Gather functions like gmt-forward in a utility-file
# decide if xp, yp etc. should be arguments so that it is
# independent of gridmap (must include optional parameters)
def gmt_forward(gmap, lon, lat):
    """Use GMT mapproject for forward projection"""

    ylon = gmap.ylon

    projection = '-JS%s/90.0/1c -C' % str(ylon)
    extent = '-R%g/%g/%g/%g' % (ylon-10, ylon+10, 60, 90)

    a, b = gmap.ellipsoid.a, gmap.ellipsoid.b
    ellipsoid = "--ELLIPSOID=%s,b=%s" % (a, b)

    m = 0.5*(1 + np.sin(gmap.lat_ts*np.pi/180))
    scale = '-F --MAP_SCALE_FACTOR=%s' % m

    program = 'GMT mapproject -bo'
    command = " ".join((program, ellipsoid, projection, extent, scale))
    print(command)
    p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE,
                     stdout=subprocess.PIPE)
    # Send lon, lat to mapproject
    # p.stdin.write("{:f} {:f}\n".format(lon, lat))
    p.stdin.write("%s %s\n" % (str(lon), str(lat)))
    out, err = p.communicate()
    x, y = struct.unpack('dd', out)
    x = x / gmap.dx + gmap.xp
    y = y / gmap.dx + gmap.yp
    return x, y


# --- Forward projection

lon, lat = 2, 66   # Station M

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
