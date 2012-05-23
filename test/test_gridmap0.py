# -*- coding: utf-8 -*-

"""Unit tests for gridmap classed"""

# ----------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# ----------------------------------

import sys
from math import pi
import unittest
import numpy as np
sys.path = ['../gridmap'] + sys.path # import from developing version
import gridmap

# ------------------------------------


class test_PolarStereographic0(unittest.TestCase):
    """Test some analytic properties of the polar stereographic map"""

    xp, yp, dx, ylon = 418.25, 257.25, 10000.0, 58.0
    map0 = gridmap.PolarStereographic(xp, yp, dx, ylon)
    map1 = gridmap.PolarStereographic(xp, yp, dx, ylon, 
                                      ellipsoid=gridmap.WGS84)

    def test_scalar(self):
        """Should return a scalar for scalar input"""
        pass

    def test_vector(self):
        """Return arrays of the same shape as the input"""
    
    def test_north_pole_forward(self):
        """The coordinates of the North Pole are xp, yp"""
        lon, lat = 17.2, 90.0
        x0, y0 = self.map0.ll2grid(lon, lat)
        self.assertEqual((x0, y0), (self.xp, self.yp))
        x1, y1 = self.map0.ll2grid(lon, lat)
        self.assertEqual((x1, y1), (self.xp, self.yp))
        
    def test_north_pole_backward(self):
        """There is no longitude at the North Pole"""
        # Should raise an exception
        lon0, lat0 = self.map0.grid2ll(self.xp, self.yp)
        lon1, lat1 = self.map1.grid2ll(self.xp, self.yp)

    def test_ylon(self):
        """lon = ylon <=> x = xp"""

        # lon = ylon => x = xp
        lon, lat = self.ylon, 72.3
        x0, y0 = self.map0.ll2grid(lon, lat)
        self.assertEqual(x0, self.xp)
        x1, y1 = self.map1.ll2grid(lon, lat)
        self.assertEqual(x1, self.xp)

        # x = xp => y = ylon
        x, y = self.xp, 222.222
        lon0, lat0 = self.map0.grid2ll(x, y)
        self.assertAlmostEqual(lon0, self.ylon, places=13)
        lon1, lat1 = self.map1.grid2ll(x, y)
        self.assertAlmostEqual(lon1, self.ylon, places=13)

    def test_inverse(self):
        """grid2ll and ll2grid are inverse"""
        lon, lat = 5.323333, 60.3925   # Bergen

        x0, y0 = self.map0.ll2grid(lon, lat)
        lon0, lat0 = self.map0.grid2ll(x0, y0)
        self.assertAlmostEqual(lon0, lon, places=14)
        self.assertEqual(lat0, lat)
        
        x1, y1 = self.map1.ll2grid(lon, lat)
        lon1, lat1 = self.map1.grid2ll(x1, y1)
        self.assertAlmostEqual(lon1, lon, places=14)
        self.assertAlmostEqual(lat1, lat, places=10)

        x, y = 200.0, 133.12345
        lon0, lat0 = self.map0.grid2ll(x, y)
        x0, y0 = self.map0.ll2grid(lon0, lat0)
        self.assertAlmostEqual(x0, x, places=12)
        self.assertAlmostEqual(y0, y, places=12)

        lon1, lat1 = self.map1.grid2ll(x, y)
        x1, y1 = self.map1.ll2grid(lon1, lat1)
        self.assertAlmostEqual(x1, x, places=9)
        self.assertAlmostEqual(y1, y, places=10)

if __name__ == '__main__':
    unittest.main()
    
        
