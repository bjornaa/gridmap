#/usr/bin/env python
# -*- coding: utf-8 -*-

"""Initiate a small demo grid"""

# Uncomment the following two lines to use the developing
import datetime
import sys
sys.path = ['..'] + sys.path
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
#Lm, Mm = 1, 1

# Name of grid
grid_name = "demo10km"
#grid_name = "demo10km_WGS84"

#
global_attributes = dict(Institution = 'Institute of Marine Research',
                         date = str(datetime.date.today()))


# -----------------------------------

# Define the grid map object, with grid size
try:
    gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, 
                   Lm=Lm, Mm=Mm, ellipsoid=ellipsoid)
except NameError:   # ellipsoid not set, using default (= sphere)
    gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, Lm, Mm)

# Create the ROMS grid file
#gridmap.create_grid(gmap, grid_name)
gridmap.create_grid(gmap, grid_name, 
                    global_attributes=global_attributes,
                    format='NETCDF4_CLASSIC')

