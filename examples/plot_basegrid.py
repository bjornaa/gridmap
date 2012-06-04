# -*- coding: utf-8 -*-

# plot the grid area using basemap

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import gridmap


# Define plotting function in basemap

def plot_basegrid(gmap, Lm, Mm):
    """Plot filled coastline for grid defined by gmap"""

    # Bør Lm og Mm være med i PolarGridMap-objektet?
    # Må nå inn som separate argument.

    if gmap.ellipsoid.name == "WGS84":
        rsphere=(gmap.ellipsoid.a, gmap.ellipsoid.b)
    else:
        rsphere = gmap.ellipsoid.a

    # Compute lon/lat of grid corners (needed by basemap=
    llcrnrlon, llcrnrlat = gmap.grid2ll(-0.5, -0.5)
    urcrnrlon, urcrnrlat = gmap.grid2ll(Lm+1.5, Mm+1.5)

    p = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                projection='stere',
                resolution='i',
                area_thresh=0.5*(gmap.dx/1000.0)**2, # km^2 
                rsphere=rsphere,
                lat_0 = 90.0,
                lon_0 = gmap.ylon,
                lat_ts = gmap.lat_c)

    p.drawcoastlines()
    p.fillcontinents(color='grey', lake_color='grey')
    # Finne på noe glupt slik at disse settes automatisk
    # ut fra lon/lat i hjørnene
    #p.drawparallels(range(50, 90, 5))
    #p.drawmeridians(range(0, 360, 15), latmax=90)
    return p

# ----------------------------

#f = Dataset('b0.nc')
#gmap = gridmap.gridmap_fromfile(f)
#Lp, Mp = len(f.dimensions['xi_rho']), len(f.dimensions['eta_rho'])
#Lm, Mm = Lp-2, Mp-2
# Grid parameters
xp    = 418.25        # x grid coordinate of north pole
yp    = 257.25        # y grid coordinate of north pole
dx    = 10000         # grid resolution (at lat_c)       [m]
ylon  = 58.0          # angle of y-axis                  [deg]

gmap = gridmap.PolarStereographic(xp, yp, dx, ylon,
                                  ellipsoid=gridmap.WGS84)

Lm, Mm = 200, 100

p = plot_basegrid(gmap, Lm, Mm)

plt.show()






