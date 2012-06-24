# -*- coding: utf-8 -*-

"""
Plot the grid area using basemap

Does not require the gridmap module
"""

# ----------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-06-23
# ----------------------------------


#import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import gridmap

# --------------------------------------
# Define the polar stereographic grid
# --------------------------------------

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at 60 )      [m]
ylon  = 58.0          # angle of y-axis                  [deg]

Lm, Mm = 200, 100     # Number of internal grid cells

ellipsoid = "WGS84"   # alternatives "sphere" or "WGS84"

# --------------------
# Other settings
# --------------------

resolution = "i" # GSHHS resolution code, "c", "l", "i", "h", or "f"
parallels = range(50, 65, 5)     # Graticule longitudes
meridians = range(-10, 20, 10)   # Graticule latitudes

# ----------------------------
# Set up the Basemap object
# ----------------------------

lat_ts = 60.0  # Latitude of true scale

if ellipsoid == "WGS84":
    a = 6378137.0
    f = 1/298.257223563
    b = (1-f)*a
    rsphere = (a, b)
else:
    rsphere = 6371000.0
    
# Basemap needs longitude, latitude of grid corners
# make a preliminary projection with origin in North Pole
p0 = Basemap(llcrnrlon=ylon, llcrnrlat=90.0,
             urcrnrlon=ylon+1, urcrnrlat=89.0,  # dummy
             projection='stere', rsphere=rsphere,
             lat_0=90.0, lon_0=ylon, lat_ts=lat_ts)

# Use full grid, including boundary cells
# extent from -0.5 to Lm+1.5 (resp. Mm+1.5)
lon_ll, lat_ll = p0((-0.5-xp)*dx, (-0.5-yp)*dx, inverse=True)
lon_ur, lat_ur = p0((Lm+1.5-xp)*dx, (Mm+1.5-yp)*dx, inverse=True)

area0 = 0.5*dx*dx/1.0e6  # Do not use islands < half grid cell

p = Basemap(llcrnrlon=lon_ll, llcrnrlat=lat_ll,
            urcrnrlon=lon_ur, urcrnrlat=lat_ur,
            projection='stere', rsphere=rsphere,
            lat_0=90.0, lon_0=ylon, lat_ts=lat_ts,
            resolution=resolution,
            area_thresh = area0)

# -------------
# Make the plot
# -------------

color='green'
p.fillcontinents(color=color, lake_color=color)

p.drawparallels(parallels)
p.drawmeridians(meridians, latmax=90)

plt.show()







