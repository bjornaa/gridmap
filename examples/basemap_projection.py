# -*- coding: utf-8 -*-

# Use basemap for a polar stereographic projection

#import numpy as np
from netCDF4 import Dataset
#import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import gridmap

# Grid parameters
xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_ts)      [m]
ylon  = 58.0          # angle of y-axis                  [deg]


lon, lat = 2, 66      # Station M
#lon, lat = 2, 90      # North Pole


# --------------------------------------------------------

# sphere case

#gmap = gridmap.PolarStereographic(xp, yp, dx, ylon,
#                        ellipsoid=gridmap.WGS84)

print "\n --- Sphere --- \n"

p = Basemap(projection='stere', rsphere=6371000,
            llcrnrlon = ylon, llcrnrlat = 90.0,
            urcrnrlon = ylon+1, urcrnrlat = 89.0,  
            lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)

x, y = p(lon, lat)
x = x / dx + xp
y = y / dx + yp

print x, y

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)
x0, y0 = gmap.ll2grid(lon, lat)
print x0, y0

print "distance [m] = ", ((x-x0)**2 + (y-y0)**2)**0.5 * dx



print "\n --- WGS84 --- \n"


gmap = gridmap.PolarStereographic(xp, yp, dx, ylon,
                                  ellipsoid=gridmap.WGS84)
x0, y0 = gmap.ll2grid(lon, lat)
print x0, y0


# WGS84 parameters
a = 6378137.0
f = 1./298.257223563
b = a*(1-f)

p = Basemap(projection='stere', rsphere=(a, b),
            llcrnrlon = ylon, llcrnrlat = 90.0,
            urcrnrlon = ylon+1, urcrnrlat = 89.0,  
            lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)

x, y = p(lon, lat)
x = x / dx + xp
y = y / dx + yp

print x, y


print "distance [m] = ", ((x-x0)**2 + (y-y0)**2)**0.5 * dx






