# -*- coding: utf-8 -*-

from __future__ import print_function

import pyproj
import gridmap


xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_ts)  [m]
ylon  = 58.0          # angle of y-axis        [deg]

lon, lat = 2, 66   # Station M
x, y = 200.0, 100.0

print("\n --- sphere ---\n")

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)

print(gmap.proj4string)



pmap = pyproj.Proj(gmap.proj4string)


x0, y0 = gmap.ll2grid(lon, lat)
x1, y1 = pmap(lon, lat)
x1, y1 = x1/gmap.dx, y1/gmap.dx

print("pyproj         : ", x1, y1)
print("gmap.ll2grid   : ", x0, y0)
print("difference [m] : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * gmap.dx)

lon0, lat0 = gmap.grid2ll(x, y)
lon1, lat1 = pmap(x*gmap.dx, y*gmap.dx, inverse=True)
print()
print("pyproj         : ", lon1, lat1)
print("gmap.grid2ll   : ", lon0, lat0)


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
