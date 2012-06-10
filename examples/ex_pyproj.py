# -*- coding: utf-8 -*-

"""Use pyproj without  gridmap to find grid coordinates of Bergen"""

import pyproj

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_ts)  [m]
ylon  = 58.0          # angle of y-axis        [deg]

lon, lat = 5.32333, 60.3925   # Bergen

# Make a dictionary of the proj parameters, WGS84 case
proj_dict = {'proj'   : 'stere',
             'ellps'  : 'WGS84',
             'lat_0'  : 90,
             'lat_ts' : 60,
             'lon_0'  : ylon,
             'x_0'    : xp*dx,
             'y_0'    : yp*dx
            }
               

pmap = pyproj.Proj(proj_dict)


x1, y1 = pmap(lon, lat)
x1, y1 = x1/dx, y1/dx

print x1, y1


lon1, lat1 = pmap(x1*dx, y1*dx, inverse=True)
print lat1

# Regner rett her, feil i test_gridmap0????
import gridmap
gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, ellipsoid=gridmap.WGS84)

print "lon, lat = ", lon, lat
x1, y1 = gmap.ll2grid(lon, lat)
print "x1, y1 = ", x1, y1
print gmap.grid2ll(x1, y1)




