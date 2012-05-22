# -*- coding: utf-8 -*-

"""Unit tests for gridmap classed"""

# ----------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# ----------------------------------

#import os
from math import pi
import unittest
import gridmap

# ------------------------------------

class test_PolarGridMap(unittest.TestCase):

    xp, yp, dx, ylon = 418.25, 257.25, 10000.0, 58.0
    map0 = gridmap.PolarGridMap(xp, yp, dx, ylon)
    map1 = gridmap.PolarGridMap(xp, yp, dx, ylon,
                                              ellipsis=gridmap.WGS84())

    def test_sphere_values(self):
        """Test some check values using spherical earth

Check values from proj4
proj +proj=stere +ellps=sphere +lat_0=90 +lon_0=ylon +lat_ts=60 \
     +x_0=xp*dx +y_0=yp*dx
"""

        lon, lat = 2, 66                   # Station M
        x0, y0 = 208.754990, 115.943832    # from proj4
        x, y = self.map0.ll2grid(lon, lat)
        self.assertAlmostEqual(x, x0, places=6) 
        self.assertAlmostEqual(y, y0, places=6) 

        lon, lat = 5.323333, 60.3925       # Bergen
        x, y = self.map0.ll2grid(lon, lat)
        x0, y0 = 168.398192, 66.753081     # from proj4
        self.assertAlmostEqual(x, x0, places=6)
        self.assertAlmostEqual(y, y0, places=6) 

    def test_WGS84_values(self):
        """Test some check values using WGS84 ellipsoid

Check values from proj4
proj +proj=stere +lat_0=90 +lon_0=ylon +lat_ts=60 +x_0=xp*dx +y_0=yp*dx
"""

        lon, lat = 2, 66                   # Station M
        x0, y0 = 207.924459, 115.383632    # from proj4
        x, y = self.map1.ll2grid(lon, lat)
        self.assertAlmostEqual(x, x0, places=6) 
        self.assertAlmostEqual(y, y0, places=6) 

        lon, lat = 5.323333, 60.3925       # Bergen
        x, y = self.map1.ll2grid(lon, lat)
        x0, y0 = 167.482134, 66.054642     # from proj4
        self.assertAlmostEqual(x, x0, places=6)
        self.assertAlmostEqual(y, y0, places=6) 
        
    def test_inverse_sphere(self):
        lon0, lat0 = 2, 66             # Station M
        x, y = self.map0.ll2grid(lon0, lat0)
        lon1, lat1 = self.map0.grid2ll(x, y)
        self.assertAlmostEqual(lon0, lon1, places=12)
        self.assertAlmostEqual(lat0, lat1, places=12)

        lon0, lat0 = 5.323333, 60.3925  # Bergen
        x, y = self.map0.ll2grid(lon0, lat0)
        lon1, lat1 = self.map0.grid2ll(x, y)
        self.assertAlmostEqual(lon0, lon1, places=12)
        self.assertAlmostEqual(lat0, lat1, places=12)

    def test_inverse_WGS84(self):
        lon0, lat0 = 2, 66             # Station M
        x, y = self.map1.ll2grid(lon0, lat0)
        lon1, lat1 = self.map1.grid2ll(x, y)
        self.assertAlmostEqual(lon0, lon1, places=12)
        self.assertAlmostEqual(lat0, lat1, places=10)

        lon0, lat0 = 5.323333, 60.3925  # Bergen
        x, y = self.map1.ll2grid(lon0, lat0)
        lon1, lat1 = self.map1.grid2ll(x, y)
        self.assertAlmostEqual(lon0, lon1, places=12)
        self.assertAlmostEqual(lat0, lat1, places=10)

    def test_pole_ll2grid(self):
        lon, lat = 23, 90      # North Pole
        x, y = self.map0.ll2grid(lon, lat)
        self.assertEqual(x, self.xp)
        self.assertEqual(y, self.yp)
        x, y = self.map1.ll2grid(lon, lat)
        self.assertEqual(x, self.xp)
        self.assertEqual(y, self.yp)

    def test_pole_grid2ll(self):
        lon, lat = self.map0.grid2ll(self.xp, self.yp)
        self.assertEqual(lat, 90.0)
        lon, lat = self.map1.grid2ll(self.xp, self.yp)
        self.assertEqual(lat, 90.0)
        
    def test_angle(self):
        
        lon, lat = self.ylon, 71.2  # at "vertical" latitude
        x, y = self.map0.ll2grid(lon, lat)
        angle = self.map0.angle(x, y)
        self.assertEqual(angle, 0.0)
        x, y = self.map1.ll2grid(lon, lat)
        angle = self.map1.angle(x, y)
        self.assertEqual(angle, 0.0)
        
        lon, lat = 2, 66   # Station M
        x, y = self.map0.ll2grid(lon, lat)
        angle = self.map0.angle(x, y)
        self.assertEqual(angle, (self.ylon-lon)*pi/180)
        x, y = self.map1.ll2grid(lon, lat)
        angle = self.map1.angle(x, y)
        self.assertEqual(angle, (self.ylon-lon)*pi/180)
        
    def test_scale(self):
        """map_scale should give 1 at 60 degrees north"""
        lon, lat = 5, 60     # at true latitude
        x, y = self.map0.ll2grid(lon, lat)
        scale = self.map0.map_scale(x, y)
        self.assertAlmostEqual(scale, 1.0, places=14)
        x, y = self.map1.ll2grid(lon, lat)
        scale = self.map1.map_scale(x, y)
        self.assertAlmostEqual(scale, 1.0, places=12)


if __name__ == '__main__':
    unittest.main()
    
        
