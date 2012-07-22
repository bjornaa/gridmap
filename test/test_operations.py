# -*- coding: utf-8 -*-

"""Unit tests for grid operations (rotation and subgrids)"""

# ----------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# ----------------------------------

import sys
import os
#from math import pi
import unittest
#import numpy as np
sys.path = ['..'] + sys.path  # import from developing version
import gridmap

class test_Projection(unittest.TestCase):

    def setUp(self):
        """Make a test projection"""
        xp, yp, dx, ylon = 418.25, 257.25, 10000, 58
        Lm, Mm = 100, 80
        self.gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, Lm, Mm)

    def test_rotate_coordinates(self):
        """Forward projection work properly"""
        gmap = gridmap.rotate(self.gmap)
        lon, lat = 2, 66  # Station M
        x0, y0 = self.gmap.ll2grid(lon, lat)
        x, y = gmap.ll2grid(lon, lat)
        self.assertAlmostEqual(x, self.gmap.Mm + 1 - y0, places=12)
        self.assertAlmostEqual(y, x0, places=12)
        
    def test_rotate4(self):
        """Rotation 4 times is the identity"""
        gmap0 = self.gmap
        gmap1 = gridmap.rotate(gmap0)
        gmap2 = gridmap.rotate(gmap1)
        gmap3 = gridmap.rotate(gmap2)
        gmap4 = gridmap.rotate(gmap3)

        self.assertAlmostEqual(gmap4.xp, gmap0.xp, places=12)
        self.assertAlmostEqual(gmap4.yp, gmap0.yp, places=12)
        self.assertAlmostEqual(gmap4.ylon, gmap0.ylon, places=12)
        
    def test_subgrid_too_large(self):
        """subgrid should fail when the subgrid is too large"""
        i0, j0 = 10, 20
        Lm = self.gmap.Lm - i0      # Maximum allowed
        Mm = self.gmap.Mm - j0      # Maximum allowed
        gridmap.subgrid(self.gmap, i0, j0, Lm, Mm) # should work
        Lm = Lm + 1   # Just too large
        self.assertRaises(ValueError, gridmap.subgrid,
                          self.gmap, i0, j0, Lm, Mm)

    def test_subgrid_coordinates(self):
        """Forward projection work properly"""
        i0, j0 = 10, 20
        Lm = self.gmap.Lm - i0 - 10
        Mm = self.gmap.Mm - j0 - 10
        gmap = gridmap.subgrid(self.gmap, i0, j0, Lm, Mm)
        lon, lat = 2, 66  # Station M
        x0, y0 = self.gmap.ll2grid(lon, lat)
        x, y = gmap.ll2grid(lon, lat)
        self.assertAlmostEqual(x, x0 - i0, places=12)
        self.assertAlmostEqual(y, y0 - j0, places=12)

# -----------------------------------

class test_GridFile(unittest.TestCase):

    def setUp(self):
        """Make a small grid file for testing"""
        xp, yp, dx, ylon = 418.25, 257.25, 10000, 58
        Lm, Mm = 30, 20
        self.grid_name = "test10km"
        self.file_name = self.grid_name + "_grid.nc"
        gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, Lm, Mm)
        gridmap.create_grid(gmap, self.grid_name, self.file_name)

    def tearDown(self):
        """Remove the grid file"""
        #os.remove(self.file_name)


    def test_rotate(self):
        pass

    def test_subgridfile(self):
        i0, j0 = 10, 5
        Lm, Mm = 10, 8
        gridmap.subgridfile(self.file_name, 'a.nc', i0, j0, Lm, Mm)
        # en test: samme antall linjer i cdl-fil (ncdump -h)v
        # test: lon m.m. er den samme

        

# ----------------------

if __name__ == '__main__':
    unittest.main()
