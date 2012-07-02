import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import gridmap

R = 6371000.0
rad = np.pi/180.0
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# --------------------

def dist(lon0, lat0, lon1, lat1):
    """Haversine formula for distances on spherical earth"""
    phi0    = lat0*rad
    phi1    = lat1*rad
    dphi    = phi1 - phi0
    dlambda = (lon1 - lon0) * rad
    a = sin(0.5*dphi)**2 + cos(phi0)*cos(phi1)*sin(0.5*dlambda)**2
    return 2 * R * np.arctan2(sqrt(a), sqrt(1-a))

# ----------------

f = Dataset('demo10km_grd.nc')
gmap = gridmap.fromfile(f)

print "xp, yp, dx, ylon = ", gmap.xp, gmap.yp, gmap.dx, gmap.ylon
print "Lm, Mm = ", gmap.Lm, gmap.Mm

pn = f.variables['pn'][:,:]
pm = f.variables['pm'][:,:]
dndx = f.variables['dndx'][:,:]
dmde = f.variables['dmde'][:,:]

lon   = f.variables['lon_rho'][:,]
lon_u = f.variables['lon_u'][:,:]
lon_v = f.variables['lon_v'][:,:]

lat   = f.variables['lat_rho'][:,]
lat_u = f.variables['lat_u'][:,:]
lat_v = f.variables['lat_v'][:,:]

# Choose a grid cell
i, j = 4, 7

# Explicit differencing
pm0 = 1.0 / dist(lon_u[j,i+1], lat_u[j,i+1], lon_u[j,i], lat_u[j,i])
pn0 = 1.0 / dist(lon_v[j+1,i], lat_v[j+1,i], lon_v[j,i], lat_v[j,i])

pm1 = gmap.map_scale(float(i),float(j)) / gmap.dx




print "lon, lat = ", lon[j,i], lat[j,i]

print ""

print "1/pm   = ", 1/pm[j,i]
print "1/pm1  = ", 1/pm1
print "1/pn0, 1/pm0 = ", 1/pm0, 1/pn0

print ""

print "dndx = ", dndx[j,i]


print "diff = ", 0.5/pn[j,i+1]-0.5/pn[j,i-1]

R = gmap.ellipsoid.a
phi0 = gmap.lat_ts*rad

m = gmap.map_scale(i,j)

dx = gmap.dx

x = (i - gmap.xp)*gmap.dx

A = - (1/m**2)  * 1 / (R*R*(1+sin(phi0))) * x * dx * dx

print "A = ", A / dx
print "A / dndx = ", A/dndx[j,i]

