# -*- coding: utf-8 -*-

### OOPS, missing France!!! (Eurasia)
### med små verdier av Lm, Mm

"""Make a coast line for a polar stereographic grid

   The polygons are given in grid coordinates,
   The polygons are cloased (useful for fill)
   and clipped against the grid boundaries

"""

# -------------------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-05-27
# -------------------------------------------


import numpy as np
from mpl_toolkits.basemap import Basemap
import gridmap

# --------------
# Functions
# --------------

def make_coast(bmap):

    X, Y = [], []

    for polygon, ptype in zip(p.coastpolygons, p.coastpolygontypes):
        if ptype == 1: # Only use land/sea boundaries
            X = X + list(polygon[0])
            Y = Y + list(polygon[1])
            X.append(np.nan)
            Y.append(np.nan)

    # Adjust for origin, in the center of the first grid cell
    X = np.array(X)/gmap.dx - 0.5
    Y = np.array(Y)/gmap.dx - 0.5

    return X, Y

def save_coast(filename, X, Y):
    fmt = '%9.4f %9.4f\n'
    f = open(filename, 'w')
    for xy in zip(X, Y):
        f.write(fmt % xy)
    f.close()

if __name__ == '__main__':

    # -------------------
    # User settings 
    # -------------------

    xp, yp, dx, ylon = 418.25, 257.25, 10000, 58.0
    Lm, Mm = 125, 100
    #Lm, Mm = 880, 660

    # Set GSHHS resolution, 'c', 'l', 'i', 'h', 'f'
    # (crude, low, intermediate, high, or full
    
    resolution = 'i'  

    gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)

    p = gmap.basemap(Lm, Mm, resolution, area_thresh=None)

    X, Y = make_coast(p)

    save_coast('a.dat', X, Y)

    
    

    






