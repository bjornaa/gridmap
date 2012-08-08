# -*- coding: utf-8 -*-

"""
Reverse engingeer a polar stereographic grid to find the grid mapping

Assumes spherical earth, R = 6361 km
and true scale at 60 degrees north

"""

import sys
import numpy as np
from netCDF4 import Dataset
import gridmap


rad = np.pi / 180.0
deg = 180.0 / np.pi


def find_parameters(f):
    """
    Find the polar stereographic parameters from a ROMS netcdf file

    Input:
      f: a netCDF4.Dataset from a grid file

    Returns a PolarStereographic instance

    The variables lon_rho, lat_rho, angle, and pm are required
    in the file.

    Return some nonsence if the grid is not polar stereographic
    with center at North pole, spherical earth radius 6371 km
    and true length at 60 degrees. The function is_polstereo can be
    used to detect this.
    """
     
    lat_ts = 60.0  # Latitude of true scale

    # The number of grid cells
    Lm = len(f.dimensions['xi_rho'])-2
    Mm = len(f.dimensions['eta_rho'])-2
    #print "Lm, Mm =", Lm, Mm

    try:
        lon0 = f.variables['lon_rho'][0,0]
        lat0 = f.variables['lat_rho'][0,0]
        angle0 = f.variables['angle'][0,0]
        pm0 = f.variables['pm'][0,0]
    except KeyError as e:
        v = "Missing essential variable %s in %s" % (e, ROMSfile)
        raise KeyError(v)
        sys.exit(1)

    ylon = lon0 + angle0*deg
    ylon = round(ylon, 2)
    dx = 1.0/pm0 * (1 + np.sin(lat_ts*rad))/(1 + np.sin(lat0*rad))
    dx = round(dx, 2)

    # Make a grid mapping with origin at north pole
    gmap0 = gridmap.PolarStereographic(0.0, 0.0, dx, ylon, 1, 1)
    x, y = gmap0.ll2grid(lon0, lat0)
    xp = round(-x, 3)
    yp = round(-y, 3)

    return gridmap.PolarStereographic(xp, yp, dx, ylon, Lm, Mm)

# --------------------------------------------

def is_polstereo(f, gmap, i, j):
    """
    Control the grid mapping of a file

    Input:       
      f : a netCDF4.Dataset from a grid file
      gmap : a PolarSterographic grid mapping
      i, j : indices of grid cell to test

    Output:
      True/False depending on the success
    """
    good = True   # no wrong tests yet

    # Read some variables at the grid cell
    lon0 = f.variables['lon_rho'][j, i]
    lat0 = f.variables['lat_rho'][j, i]
    angle0 = f.variables['angle'][j, i]

    # Improve tests, find reasonable tolerances
    tol = 1.0e-6
    lon, lat = gmap.grid2ll(i, j)
    good = good and (abs(lon-lon0) < tol) and (abs(lat-lat0) < tol)

    good = good and (abs(gmap.angle(i,j)-angle0) < tol)

    return good

# ----------------------------------------------

if __name__ == '__main__':

    import sys
    
    #usage = "Usage: python find_parameters ROMSfile"
    usage = "Usage: python %s ROMSfile" % sys.argv[0]

    try:
        ROMSfile = sys.argv[1]
    except IndexError:
        print usage
        sys.exit(1)

    try:
        f = Dataset(ROMSfile)
    except RuntimeError as e:
        print "%s : %s" % (ROMSfile, e)
        sys.exit(1)
 
    # Find the parameteres
    try:
        gmap = find_parameters(f)
    except KeyError as e:
        print "KeyError:", e
        sys.exit(1)
        
    i = gmap.Lm // 2
    j = gmap.Mm // 3

    good = is_polstereo(f, gmap, j, i)
    
    if not good:
        print "%s : not polar stereographic grid mapping" % ROMSfile
        sys.exit(1)

    print "%s : polar stereographic grid mapping" % ROMSfile
    print "xp   = ", gmap.xp
    print "yp   = ", gmap.yp
    print "dx   = ", gmap.dx
    print "ylon = ", gmap.ylon
    print "Lm   = ", gmap.Lm
    print "Mm   = ", gmap.Mm



        
        

    
