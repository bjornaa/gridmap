# -*- coding: utf-8 -*-


"""Make a coast line for a polar stereographic grid

   The polygons are given in grid coordinates,
   The polygons are cloased (useful for fill)
   and clipped against the grid boundaries

"""

# Note: Bug in 1.0.3  european continent may be missing
# Fixed in 1.0.4

# -------------------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-05-27
# -------------------------------------------

import numpy as np
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import gridmap

# --------------
# Functions
# --------------

def get_coast(bmap):
    """Use basemap to extract the coast line

    Input:
        bmap: Basemap instance.
    Output:
        X, Y: 1D arrays of x, y coordinates of coast line.
              Segments are separated by Nan's
    """
    
    X, Y = [], []

    for polygon, ptype in zip(p.coastpolygons, p.coastpolygontypes):
        if ptype == 1: # Only use land/sea boundaries
            X = X + list(polygon[0])
            Y = Y + list(polygon[1])
            X.append(np.nan)
            Y.append(np.nan)

    return np.array(X), np.array(Y)

def save_coast(filename, X, Y):
    """Save a coast line to a text file"""
    fmt = '%9.4f %9.4f\n'
    f = open(filename, 'w')
    for xy in zip(X, Y):
        f.write(fmt % xy)
    f.close()

if __name__ == '__main__':

    # -------------------
    # User settings 
    # -------------------

    gridfile  = "demo10km_grid.nc"
    coastfile = "demo10km_coast.dat"

    # Set GSHHS resolution, 'c', 'l', 'i', 'h', 'f'
    # (crude, low, intermediate, high, or full
    resolution = 'i'  

    # Define gridmap instance from the grid file
    f = Dataset(gridfile)
    gmap = gridmap.fromfile(f)

    # Define a Basemap instance
    p = gmap.basemap(resolution)

    # Make basemap extract the coast line
    X, Y = get_coast(p)
    # Scale to grid coordinates
    # NB. What should be origin in gridmap.basemap???
    X = X/gmap.dx
    Y = Y/gmap.dx

    save_coast(coastfile, X, Y)

    
    

    






