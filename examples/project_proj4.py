# -*- coding: utf-8 -*-

#import sys
import subprocess
import struct
from math import pi
from gridmap import PolarStereographic


# -----------------------------------------------------

def proj(projstring, lon, lat):
    """Use proj4 for forward projection"""

    # Use -o for binary outpur => no errors due to format
    command = 'proj -o ' + projstring

    # Set up the proj process
    #if verbose: print command
    print command
    p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

    # Send lon, lat to proj
    p.stdin.write("%s %s\n" % (str(lon), str(lat)))

    # Get the output 
    out, err = p.communicate()
    x, y = struct.unpack('2d', out)
    return x, y

def invproj(projstring, x, y):
    """Use proj4 for inverse projection"""

    command = 'invproj -o ' + projstring

    #if verbose: print command
    p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

    p.stdin.write("%s %s\n" % (str(x), str(y)))

    out, err = p.communicate()
    lon, lat = struct.unpack('2d', out)
    # Convert from radians to degrees
    lon = lon * 180.0 / pi
    lat = lat * 180.0 / pi

    return lon, lat

# ----------------------------------------------------

if __name__ == '__main__':

    xp    = 418.25        # x grid coordinate of north pole
    yp    = 257.25        # y grid coordinate of north pole
    dx    = 10000         # grid resolution (at lat_ts)  [m]
    ylon  = 58.0          # angle of y-axis        [deg]

    lon, lat = 2, 66   # Station M
    x, y = 0.0, 0.0

    #verbose = True

    print "\n --- sphere ---\n"

    gmap = PolarStereographic(xp, yp, dx, ylon)

    x0, y0 = gmap.ll2grid(lon, lat)
    x1, y1 = proj(gmap.proj4string, lon, lat)
    x1, y1 = x1/gmap.dx, y1/gmap.dx

    print "proj           : ", x1, y1
    print "gmap.ll2grid   : ", x0, y0
    print "difference [m] : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * gmap.dx

    lon0, lat0 = gmap.grid2ll(x, y)
    lon1, lat1 = invproj(gmap.proj4string, x*gmap.dx, y*gmap.dx)
    print
    print "invproj        : ", lon1, lat1
    print "gmap.grid2ll   : ", lon0, lat0


    print "\n --- WGS84 ---\n"

    gmap = PolarStereographic(xp, yp, dx, ylon, ellipsoid='WGS84')

    x0, y0 = gmap.ll2grid(lon, lat)
    x1, y1 = proj(gmap.proj4string, lon, lat)
    x1, y1 = x1/gmap.dx, y1/gmap.dx

    print "proj           : ", x1, y1
    print "gmap.ll2grid   : ", x0, y0
    print "difference [m] : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * gmap.dx

    gmap = PolarStereographic(xp, yp, dx, ylon, ellipsoid='WGS84')
    lon0, lat0 = gmap.grid2ll(x, y)
    lon1, lat1 = invproj(gmap.proj4string, x*gmap.dx, y*gmap.dx)
    print
    print "invproj        : ", lon1, lat1
    print "gmap.grid2ll   : ", lon0, lat0








