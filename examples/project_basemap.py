# -*- coding: utf-8 -*-

"""Use basemap for a polar stereographic projection"""

# This example uses gridmap only for comparison

# --------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-06-29
# --------------------------------

from mpl_toolkits.basemap import Basemap
import gridmap

# Grid parameters
xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_ts)      [m]
ylon  = 58.0          # angle of y-axis                  [deg]

# Position to project
lon, lat = 2, 66      # Station M

# ------------------------
# Sphere case
# ------------------------

print "\n --- Sphere --- \n"

# Set up projection in basemap, origin at North Pole
# with dummy values for urcrnrlon and urcrnrlat
bmap = Basemap(projection='stere', rsphere=6371000,
               llcrnrlon = ylon, llcrnrlat = 90.0,
               urcrnrlon = ylon+1, urcrnrlat = 89.0,  
               lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)

x, y = bmap(lon, lat)
# Rescale to grid coordinates
x = x / dx + xp
y = y / dx + yp

# Set up projection in gridmap
gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)
x0, y0 = gmap.ll2grid(lon, lat)

print "basemap: ", x, y
print "gridmap: ", x0, y0
print "distance [m] = ", ((x-x0)**2 + (y-y0)**2)**0.5 * dx

# ------------------------
# WGS84 case
# ------------------------

print "\n --- WGS84 --- \n"

# WGS84 parameters for basemap
a = 6378137.0
f = 1./298.257223563
b = a*(1-f)

# Projection in basemap
bmap = Basemap(projection='stere', rsphere=(a, b),
            llcrnrlon = ylon, llcrnrlat = 90.0,
            urcrnrlon = ylon+1, urcrnrlat = 89.0,  
            lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)

x, y = bmap(lon, lat)
x = x / dx + xp
y = y / dx + yp

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, ellipsoid='WGS84')
x0, y0 = gmap.ll2grid(lon, lat)

print "basemap: ", x, y
print "gridmap: ", x0, y0
print "distance [m] = ", ((x-x0)**2 + (y-y0)**2)**0.5 * dx

