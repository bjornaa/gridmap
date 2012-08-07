# -*- encoding: utf-8 -*-

"""
Reverse engingeer a polar stereographic grid to find the grid mapping

Assumes spherical earth, R = 6361 km
and true scale at 60 degrees north

"""

import numpy as np
from netCDF4 import Dataset
import gridmap

gridfile = "demo10km_grid.nc"

rad = np.pi / 180.0
deg = 180.0 / np.pi

#xp0 = 0.0
#yp0 = 0.0
#Lm0 = 1
#Mm0 = 1

f = Dataset(gridfile)

# The number of grid cells
Lm = len(f.dimensions['xi_rho'])-2
Mm = len(f.dimensions['eta_rho'])-2
print "Lm, Mm =", Lm, Mm


# lon_rho
lon0 = f.variables['lon_rho'][0,0]
angle0 = f.variables['angle'][0,0]
ylon = lon0 + angle0*deg
print "ylon = ", ylon

# dx
pm0 = f.variables['pm'][0,0]
lat0 = f.variables['lat_rho'][0,0]
dx0 = 1.0/pm0
dx = dx0 * (1 + np.sin(60*rad))/(1 + np.sin(lat0*rad))
# Burde kanskje runde av lite grann
print "dx = ", dx





# Make a grid with origin at north pole
gmap0 = gridmap.PolarStereographic(0.0, 0.0, dx, ylon, 1, 1)


#lon = f0.variables['lon_rho'][:,:]
#lat = f0.variables['lat_rho'][:,:]
#angle0 = f0.variables['angle'][0,0]


# gmap0 coordinates of origo
x, y = gmap0.ll2grid(lon0, lat0)

xp, yp = -x, -y

print "xp, yp = ", xp, yp

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, Lm, Mm)

# --------
# Control
# --------

# OK, origin in place
i, j = 0, 0
lo, la = gmap.grid2ll(j, i)
print "--- i,j = ", i, j
print "reconstructed grid : ", lo, la
print "in file            : ", lon[j,i], lat[j,i]

i, j = 70, 50  # OK, this fits
lo, la = gmap.grid2ll(i, j)
print "--- i,j = ", i, j
print "reconstructed grid : ", lo, la
print "in file            : ", lon[j,i], lat[j,i]

# This is OK, meaning the value of ylon is OK
print "--- angle at origin [degrees]"
print "reconstructured grid : ", ylon - lon[0,0]
print "from file            : ", angle0 * 180.0/np.pi

# --------------------------
# Generate a new grid file
# --------------------------

# Commented out, enough to do once
#grid_name = "norkyst_800m"
#gridmap.create_grid(gmap, grid_name)

f = Dataset("norkyst_800m_grid.nc", 'a')

# -----------
# More tests
# -----------

i,j = 70, 50

# Should be almost equal
print "pm    : ", f0.variables['pm'][j,i], f.variables['pm'][j,i]

# Should be equal
print "lat_u : ", f0.variables['lat_u'][j,i], f.variables['lat_u'][j,i]

# -------------------------------
# Transfer topography and masks to new grid file
# -------------------------------

H = f0.variables['h'][:,:]
f.variables['hraw'][0,:,:] = H
f.variables['h'][:,:] = H

f.variables['mask_rho'][:,:] = f0.variables['mask_rho'][:,:]
f.variables['mask_u'][:,:] = f0.variables['mask_u'][:,:]
f.variables['mask_v'][:,:] = f0.variables['mask_v'][:,:]
f.variables['mask_psi'][:,:] = f0.variables['mask_psi'][:,:]

f.close()

