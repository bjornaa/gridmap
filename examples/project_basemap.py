#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Use basemap for a polar stereographic projection
# without gridmap

# ----------------------------------
# project_basemap.py
#
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-06-04
# ----------------------------------


from mpl_toolkits.basemap import Basemap

# -------------------
# User settings
# -------------------

# Grid parameters
xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_c)       [m]
ylon  = 58.0          # angle of y-axis                  [°N]

lon_ts = 60.0         # Latitude of true scale           [°N]

lonlat_file = "lonlat.dat"  # File with longitude-latitude data


# -----------------
# Read the data 
# -----------------

f = open(lonlat_file)

lonlats = []
for line in f:
    w = line.split()
    lonlats.append((float(w[0]), float(w[1])))
    
f.close()

# --------------
# Sphere case
# --------------

print " --- Sphere --- "

# Setup the Basemap projection
# urcrnlon, urcrnrlat is arbitrary
p = Basemap(projection='stere', rsphere=6371000,
            llcrnrlon = ylon, llcrnrlat = 90.0,
            urcrnrlon = ylon+1, urcrnrlat = 89.0,  
            lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)

# Compute grid coordinates for the lon, lat pairs
for lon, lat in lonlats:
    x, y = p(lon, lat)
    x = x / dx + xp
    y = y / dx + yp
    print "%9.4f %9.4f" % (x, y)

# --------
# WGS84
# --------

print " --- WGS84 --- "

# WGS84 parameters
a = 6378137.0           # Major semiaxis [m]
f = 1./298.257223563    # Flatening
b = a*(1-f)             # Minor semiaxis [m]

p = Basemap(projection='stere', rsphere=(a, b),
            llcrnrlon = ylon, llcrnrlat = 90.0,
            urcrnrlon = ylon+1, urcrnrlat = 89.0,  
            lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)

for lon, lat in lonlats:
    x, y = p(lon, lat)
    x = x / dx + xp
    y = y / dx + yp
    print "%9.4f %9.4f" % (x, y)







