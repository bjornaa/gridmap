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

