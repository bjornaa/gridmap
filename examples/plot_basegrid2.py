# -*- coding: utf-8 -*-

"""
Plot the grid area using basemap and gridmap
"""

# ----------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-06-24
# ----------------------------------

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from gridmap import PolarStereographic

# --------------------------------------
# Define the polar stereographic grid
# --------------------------------------

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at 60 )      [m]
ylon  = 58.0          # angle of y-axis                  [deg]

Lm, Mm = 200, 100     # Number of internal grid cells
Lm, Mm = 40, 30

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

gmap = PolarStereographic(xp, yp, dx, ylon, Lm, Mm, ellipsoid)
#gmap = PolarStereographic(xp+0.5, yp+0.5, dx, ylon, Lm+1, Mm+1, ellipsoid)

bmap = gmap.basemap(resolution=resolution)
# Note: bmap gives same projection as gmap
# side effect: the plot differs from plot_basegrid.py
#              by showing only half the boundary cells
#              i.e., extent from 0.0 to Lm+1 resp. Mm+1
# TODO: give basemap method an optional argument,
#       gmap.basemap(boundary=full/half/none)
# Alternative, use
# PolarStereographic(xp+0.5*dx, yp+0.5*dx, ..., Lm+1.0, Mm+1.0)

# -------------
# Make the plot
# -------------

color='green'
bmap.fillcontinents(color=color, lake_color=color)

bmap.drawparallels(parallels)
bmap.drawmeridians(meridians, latmax=90)

plt.show()







