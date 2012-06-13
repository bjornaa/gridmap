from os import system
from netCDF4 import Dataset
import gridmap

#f = Dataset('b0.nc')
#gmap = gridmap.fromfile(f)
#L = len(f.dimensions['xi_rho'])
#M = len(f.dimensions['eta_rho'])

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_c)       [m]
ylon  = 58.0          # angle of y-axis                  [deg]

Lm, Mm = 125, 100
L, M = Lm+1, Mm+1

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)

lon0, lat0 = gmap.grid2ll(0.0, 0.0)
lon1, lat1 = gmap.grid2ll(L, M)

projection = '-JS%s/90.0/22c' % str(ylon)
extent = '-R%g/%g/%g/%gr' % (lon0, lat0, lon1, lat1)
boundary = "-B30g15/20g5"
landcolor = '-Ggreen'

command = "pscoast %s %s %s %s -Wthinnest -A200  > a.ps" % \
            (projection, extent, boundary, landcolor)
print command

system(command)

system('gv a.ps')
