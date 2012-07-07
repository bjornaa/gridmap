# -*- coding: utf-8 -*-

"""Unit tests for gridmap classed"""

# ----------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# ----------------------------------

import sys
#from math import pi
import os
import unittest
import numpy as np
from netCDF4 import Dataset
sys.path = ['..'] + sys.path # import from developing version
import gridmap

rad = np.pi / 180.0
sin, cos, sqrt = np.sin, np.cos, np.sqrt
R = 6371000.0

def dist(lon0, lat0, lon1, lat1):
    """Haversine formula for distances on spherical earth"""
    phi0    = lat0*rad
    phi1    = lat1*rad
    dphi    = phi1 - phi0
    dlambda = (lon1 - lon0) * rad
    a = sin(0.5*dphi)**2 + cos(phi0)*cos(phi1)*sin(0.5*dlambda)**2
    return 2 * R * np.arctan2(sqrt(a), sqrt(1-a))


class test_GridGeneration_sphere(unittest.TestCase):
    
    def setUp(self):
        """Make a small grid file for testing"""
        xp, yp, dx, ylon = 418.25, 257.25, 10000, 58
        Lm, Mm = 10, 8
        self.grid_name = "test10km"
        self.file_name = self.grid_name + "_grid.nc"
        self.gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, Lm, Mm)
        gridmap.create_grid(self.gmap, self.grid_name, self.file_name)
        #print "setUp ferdig"



#    def tearDown(self):
#        """Remove the grid file"""
#        os.remove(self.file_name)


    def test_angle(self):
        """Angle variable is OK"""

        i, j = 4, 3
        f = Dataset(self.file_name)
        # Compute angle variable by differencing
        # like Kate Hedstrom's  gridpak
        pass

    def test_pm(self):
        """pm and pn are OK"""
        i, j = 4, 3
        f = Dataset(self.file_name)
        gmap = gridmap.fromfile(f)
        pm = f.variables['pm'][:,:]
        pn = f.variables['pn'][:,:]
        lon_u = f.variables['lon_u'][:,:]
        lat_u = f.variables['lat_u'][:,:]
        lon_v = f.variables['lon_v'][:,:]
        lat_v = f.variables['lat_v'][:,:]
        f.close()

        # Difference formulation, using spherical distance
        pm0 = 1.0 / dist(lon_u[j,i+1], lat_u[j,i+1], lon_u[j,i], lat_u[j,i])
        pn0 = 1.0 / dist(lon_v[j+1,i], lat_v[j+1,i], lon_v[j,i], lat_v[j,i])

        # Analytical from map_scale
        pm1 = gmap.map_scale(float(i),float(j)) / gmap.dx

        self.assertAlmostEqual(pm[j,i], pm0, places=6)
        self.assertAlmostEqual(pn[j,i], pn0, places=6)
        # Do not expect equality if variables saved as float32
        self.assertEqual(pm[j,i], pm1)
        self.assertEqual(pn[j,i], pm1)
    
    def test_dmde(self):
        """dmde and dndx are OK"""


        # Hent kode fra examples/aa.py
        #   som kan slettes etterpå
        i, j = 4, 3
        f = Dataset(self.file_name)
        gmap = gridmap.fromfile(f)
        pm = f.variables['pm'][:,:]
        pn = f.variables['pn'][:,:]
        lon_u = f.variables['lon_u'][:,:]
        lat_u = f.variables['lat_u'][:,:]
        lon_v = f.variables['lon_v'][:,:]
        lat_v = f.variables['lat_v'][:,:]
        f.close()


    def test_fromfile(self):
        """Testing that fromfile recreates the grid map object"""
        f = Dataset(self.file_name)
        gmap1 = gridmap.fromfile(f)
        gmap0 =self.gmap
        for att in ['xp', 'yp', 'dx', 'ylon', 'Lm','Mm', 'lat_ts']:
            self.assertEqual(getattr(gmap0, att), getattr(gmap1, att))
        self.assertEqual(gmap1.ellipsoid.a, gmap0.ellipsoid.a)
        self.assertEqual(gmap1.ellipsoid.invf, gmap0.ellipsoid.invf)
        
if __name__ == '__main__':
    unittest.main()

