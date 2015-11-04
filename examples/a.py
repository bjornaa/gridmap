# -*- coding: utf-8 -*-

from __future__ import print_function

from os import system
import subprocess
import struct

from netCDF4 import Dataset
import gridmap

#f = Dataset('b0.nc')
#gmap = gridmap.fromfile(f)
#L = len(f.dimensions['xi_rho'])
#M = len(f.dimensions['eta_rho'])

xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_ts)      [m]
ylon  = 58.0          # angle of y-axis                  [deg]

Lm, Mm = 125, 100
L, M = Lm+1, Mm+1

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)

lon0, lat0 = gmap.grid2ll(0.0, 0.0)
lon1, lat1 = gmap.grid2ll(L, M)

projection = '-JS%s/90.0/%sc' % (str(ylon), str(L))
extent = '-R%g/%g/%g/%gr' % (lon0, lat0, lon1, lat1)
# scale = '--PROJ_SCALE_FACTOR=2.2'
scale = '-F'
ellipsoid = "--PROJ_ELLIPSOID=%s" % "6371000,b=6371000"

command0 = 'gmt mapproject -bo'
# gmtstring = " ".join((ellipsoid, projection, offset, extent))
command = " ".join((command0,  projection, extent, scale, ellipsoid))

print(command)
p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE,
                     stdout=subprocess.PIPE)

# Send lon, lat to mapproject
# p.stdin.write("{:f} {:f}\n".format(lon, lat))
# p.stdin.write("%s %s\n" % (str(lon0), str(lat0)))
p.stdin.write("%s %s\n" % (str(lon1), str(lat1)))
# p.stdin.write("%s %s\n" % ('0', '90'))


# Get the output
out, err = p.communicate()
x, y = struct.unpack('dd', out)

# x = x / dx
# y = y / dx
# x = L*x
# y = L*y

print(x, y)
