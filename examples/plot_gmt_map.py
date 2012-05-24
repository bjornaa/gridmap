from os import system
from netCDF4 import Dataset
import gridmap

f = Dataset('b0.nc')

v = f.variables['arctic10_WGS84'] # grid mapping definition
if v.ellipsoid == 'WGS84':
    ellipsoid = gridmap.WGS84()
else:
    ellipsoid = gridmap.Sphere(v.earth_radius)
xp = v.false_easting
yp = v.false_northing
ylon = v.straight_vertical_longitude_from_pole
dx = v.dx

gmap = gridmap.PolarGridMap(xp,yp,dx,ylon,ellipsoid=ellipsoid)

L = len(f.dimensions['xi_rho'])
M = len(f.dimensions['eta_rho'])

lon0, lat0 = gmap.grid2ll(0.0, 0.0)
lon1, lat1 = gmap.grid2ll(L, M)

projection = '-JS%s/90.0/22c' % str(ylon)
extent = '-R%g/%g/%g/%gr' % (lon0, lat0, lon1, lat1)
boundary = "-B30g15/20g5"
landcolor = '-Ggreen'

command = "/usr/lib/gmt/bin/pscoast %s %s %s %s -Wthinnest -A200  > a.ps" % \
            (projection, extent, boundary, landcolor)
print command

system(command)

system('gv a.ps')
