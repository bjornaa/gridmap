# -*- coding: utf-8 -*-

"""
Example using pyproj without gridmap for polar stereogaphic grid projection

The example requires pyproj or basemap
pyproj: https://code.google.com/p/pyproj/
basemap: http://matplotlib.github.com/basemap/

"""

# ----------------------------------
# ex_pyproj.py
#
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# Bergen, Norway
# 2012-06-10
# -----------------------------------

from __future__ import print_function

#import pyproj                                 # separate pyproj
import mpl_toolkits.basemap.pyproj as pyproj  # pyproj from basemap

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_ts)  [m]
ylon  = 58.0          # angle of y-axis        [deg]

lon, lat = 5.32333, 60.3925   # Bergen

# -----------------
# sphere case
# -----------------

# Make a dictionary of the proj parameters
proj_dict = {'proj'   : 'stere',
             'R'      : 6371000,
             'lat_0'  : 90,
             'lat_ts' : 60,
             'lon_0'  : ylon,
             'x_0'    : xp*dx,
             'y_0'    : yp*dx
            }

# Make a projection instance
pmap = pyproj.Proj(proj_dict)


x, y = pmap(lon, lat)
x, y = x/dx, y/dx
print("sphere: x, y = ", x, y)

lon1, lat1 = pmap(x*dx, y*dx, inverse=True)
print("sphere: lon, lat = ", lon1, lat1)

# -------------
# WGS84 case
# -------------

# Modify the dictionary for WGS84
proj_dict.pop('R')              # Remove radius field
proj_dict['ellps'] = 'WGS84'    # Add ellipsoid field

pmap = pyproj.Proj(proj_dict)

x, y = pmap(lon, lat)
x, y = x/dx, y/dx
print("WGS84:  x, y = ", x, y)

lon1, lat1 = pmap(x*dx, y*dx, inverse=True)
print("WGS84:  lon, lat = ", lon1, lat1)






