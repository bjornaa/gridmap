# -*- coding: utf-8 -*-

"""Make a coast line for a polar stereographic grid

   The polygons are given in grid coordinates,
   The polygons are cloased (useful for fill)
   and clipped against the grid boundaries

"""

import numpy as np
from netCDF4 import Dataset
#import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import gridmap

# -------------------
# User settings 
# -------------------

gridfile = 'b0.nc'

# Set GSHHS resolution, 'c', 'l', 'i', 'h', 'f'
# (crude, low, intermediate, high, or full
resolution = 'i'  
#area_thresh = 60  # areal threshold [km^2]

# Open grid file
f = Dataset('b0.nc')

# Map definition variable
v = f.variables['grid_mapping']

if v.ellipsoid == 'WGS84':
    ellipsoid = gridmap.WGS84()
    rsphere = (ellipsoid.a, ellipsoid.b)
else:
    ellipsoid = gridmap.Sphere(v.earth_radius)
    rshpere = ellipsoid.a

xp   = v.false_easting
yp   = v.false_northing
ylon = v.straight_vertical_longitude_from_pole
dx   = v.dx
#lat_c = v.standard_parallel
lat_c = 60

gmap = gridmap.PolarGridMap(xp,yp,dx,ylon,ellipsoid=ellipsoid)
area_thresh = 0.5*(dx*dx) / 1.e6  # unit = square kilometer

Lp, Mp = len(f.dimensions['xi_rho']), len(f.dimensions['eta_rho'])

# Compute lon/lat of grid corners (needed by basemap)
lon0, lat0 = gmap.grid2ll(-0.5, -0.5)
lon1, lat1 = gmap.grid2ll(Lp-0.5, Mp-0.5)

p = Basemap(llcrnrlon=lon0, llcrnrlat=lat0,   # lower-left corner
            urcrnrlon=lon1, urcrnrlat=lat1,   # upper-right corner
            projection='stere',           
            resolution=resolution,
            area_thresh=area_thresh,  # km^2 oppløsningsavhengig
            rsphere=rsphere,
            lat_0 = 90.0,
            lon_0 = gmap.ylon,
            lat_ts = gmap.lat_c)

X = []
Y = []

print "lager polygon"
for polygon, ptype in zip(p.coastpolygons, p.coastpolygontypes):
    if ptype == 1: # Only use land/sea boundaries
        X = X + list(polygon[0])
        Y = Y + list(polygon[1])
        X.append(np.nan)
        Y.append(np.nan)

# Adjust for origin, in the center of the first grid cell
X = np.array(X)/dx - 0.5
Y = np.array(Y)/dx - 0.5

# Liten test: ser at transform er like
#print gmap.ll2grid(5, 60)
#print [s/dx-0.5 for s in p(5, 60)]
#print "\n"
#print gmap.ll2grid(30, 90)
#print [s/dx-0.5 for s in p(30, 90)]

print "saving coast file"
f0 = open(f.gridname + '.dat', 'w')
for x, y in zip(X, Y):
    f0.write("%g %g\n" % (x, y))
f0.close()


#plt.fill(X, Y)

#plt.axis('image')

#plt.show()








