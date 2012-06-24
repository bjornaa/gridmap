# -*- coding: utf-8 -*-

#import sys
import subprocess
import struct
from math import pi
from gridmap import PolarStereographic

#verbose=True

# -----------------------------------------------------

def mapproject(gmtstring, lon, lat, verbose=False):
    """Use gmt mapproject for forward projection"""

    # Use -bo for binary outpur => no errors due to format
    command = 'mapproject -bo ' + gmtstring

    # Set up the  process
    if verbose: print command 
    p = subprocess.Popen(command, shell=True,
                     stdin=subprocess.PIPE, 
                     stdout=subprocess.PIPE)

    # Send lon, lat to mapproject
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
    x, y = 200.0, 100.0

    #verbose = True

    # ----------------------------
    print "\n --- sphere ---\n"
    # ----------------------------

    gmap = PolarStereographic(xp, yp, dx, ylon)

    projection = '-Js%s/90.0/%s/1:%s' % \
                 (str(ylon), str(gmap.lat_ts), str(100*dx))
    #ellipsoid = "--ELLIPSOID=Sphere"
    #ellipsoid = "--ELLIPSOID=WGS-84"
    ## Ser ut som ellipsoid er gitt på korrekt vis, men får feil
    ## Setting virker for WGS84 
    ellipsoid = "--ELLIPSOID=%s" % str(gmap.ellipsoid.a)
    #ellipsoid = "--ELLIPSOID=%s" % "6378137,298.257223563" # WGS84 OK
    #ellipsoid = "--ELLIPSOID=%s" % "6371000,f=0"

    extent = '-R0/1/60/61'   # Actual values are not used
    offset = '-C%s/%s' % (str(xp), str(yp))
    gmtstring = " ".join((ellipsoid, projection, offset, extent))

    x0, y0 = gmap.ll2grid(lon, lat)
    x1, y1 = mapproject(gmtstring, lon, lat, verbose=True)


    print "gmt mapproject : ", x1, y1
    print "gmap.ll2grid   : ", x0, y0
    print "difference [m] : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * gmap.dx
    print
    
    gmap = PolarStereographic(xp, yp, dx, ylon)
    lon0, lat0 = gmap.grid2ll(x, y)
    lon1, lat1 = mapproject('-I ' + gmtstring, x, y, verbose=True)
    print "mapproject -I  : ", lon1, lat1
    print "gmap.grid2ll   : ", lon0, lat0


    # ---------------------------
    print "\n --- WGS84 ---\n"
    # ---------------------------

    gmap = PolarStereographic(xp, yp, dx, ylon, ellipsoid="WGS84")

    projection = '-Js%s/90.0/%s/1:%s' % \
                 (str(ylon), str(gmap.lat_ts), str(100*dx))
    extent = '-R0/1/60/61'   # Actual values are not used
    offset = '-C%s/%s' % (str(xp), str(yp))
    gmtstring = " ".join((projection, offset, extent))


    x0, y0 = gmap.ll2grid(lon, lat)
    x1, y1 = mapproject(gmtstring, lon, lat, verbose=True)


    print "gmt mapproject : ", x1, y1
    print "gmap.ll2grid   : ", x0, y0
    print "difference [m] : ", ((x1-x0)**2 + (y1-y0)**2)**0.5 * gmap.dx
    print
    
    gmap = PolarStereographic(xp, yp, dx, ylon, ellipsoid="WGS84")
    lon0, lat0 = gmap.grid2ll(x, y)
    lon1, lat1 = mapproject('-I ' + gmtstring, x, y, verbose=True)
    print "mapproject -I  : ", lon1, lat1
    print "gmap.grid2ll   : ", lon0, lat0

