import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import gridmap

f = Dataset('demo10km_grd.nc')
gmap = gridmap.fromfile(f)

dndx = f.variables['dndx'][:,1:-1]

Lm, Mm = gmap.Lm, gmap.Mm

Lvert = 2*Lm + 5
Mvert = 2*Mm + 5
X0 = 0.5*np.arange(Lvert)-0.5
Y0 = 0.5*np.arange(Mvert)-0.5
# Make 2D arrays with grid coordonates
Xvert, Yvert = np.meshgrid(X0, Y0)
Xrho = Xvert[1::2, 1::2]
Yrho = Yvert[1::2, 1::2]

lat = f.variables['lat_rho'][:,1:-1]

print "xp, yp, dx, ylon = ", gmap.xp, gmap.yp, gmap.dx, gmap.ylon

rad = np.pi / 180.0


pn = gmap.map_scale(Xrho, Yrho) / gmap.dx

dndx2 = 0.5/pn[:,2:] - 0.5/pn[:,:-2]

X = Xrho[:,1:-1]
B = (X - gmap.xp) * gmap.dx * (1 + np.sin(lat*rad))**2
B = B / (1+np.sin(gmap.lat_ts*rad))**3

# B / dndx2 er på det nærmeste konstant
# men forksjell er systematisk

