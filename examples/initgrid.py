#/usr/bin/env python
# -*- coding: utf-8 -*-

"""Initiate a small demo grid"""

import gridmap

# ----------------
# User settings
# ----------------

# Polar stereographic map projection
xp, yp, dx, ylon = 418.25, 257.25, 10000, 58
# Uncomment the following line for WGS84 ellipsoid
#ellipsoid = 'WGS84'

# Number of internal grid cells
# deliberately small numbers for the demo
Lm, Mm = 100, 75

# Name of grid
grid_name = "demo10km"

# Name of grid file
file_name = grid_name + "_grd.nc"

# -----------------------------------

# Define the grid map object, with grid size
try:
    gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, 
                   Lm=Lm, Mm=Mm, ellipsoid=ellipsoid)
except NameError:   # ellipsoid not set, using default (= sphere)
    gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, Lm, Mm)

# Create the ROMS grid file
gridmap.create_grid(gmap, grid_name, file_name)

