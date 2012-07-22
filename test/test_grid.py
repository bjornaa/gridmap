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
        #gridmap.create_grid(self.gmap, self.grid_name, self.file_name)
        gridmap.create_grid(self.gmap, self.grid_name)
        #print "setUp ferdig"



    def tearDown(self):
        """Remove the grid file"""
        os.remove(self.file_name)

    # ------------------------------------------------

    def test_metric(self):
        """Check the (geo)metric variables

        They should agree with finite difference approximation
        as done in Kate Hedstrom's gridpak.
        They should agree with analytical formulation for
        the polar stereographic map projection

        Presently the test only applies to the spherical case
        """

        # Take an "arbitrary" grid cell
        i, j = 4, 3

        # Read the grid file
        f = Dataset(self.file_name)       
        gmap = gridmap.fromfile(f)
        angle = f.variables['angle'][:,:]
        pm = f.variables['pm'][:,:]
        pn = f.variables['pn'][:,:]
        dmde = f.variables['dmde'][:,:]
        dndx = f.variables['dndx'][:,:]
        lon_rho = f.variables['lon_rho'][:,:]
        lon_u = f.variables['lon_u'][:,:]
        lat_u = f.variables['lat_u'][:,:]
        lon_v = f.variables['lon_v'][:,:]
        lat_v = f.variables['lat_v'][:,:]
        f.close()

        # Compute angle variable by differencing

        # U-points
        a1 = lat_u[j,i+1] - lat_u[j,i]
        a2 = lon_u[j,i+1] - lon_u[j,i]
        #a2 = np.where(a2 < -180., a2+360., a2)
        #a2 = np.where(a2 > +180., a2-360., a2)
        a2 = a2 * cos(0.5*(lat_u[j,i+1] + lat_u[j,i])*rad)
        angle_u = np.arctan2(a1,a2)

        # V-points
        a2 = lat_v[j,i] - lat_v[j+1, i]
        a1 = lon_v[j,i] - lon_v[j+1, i]
        #a1 = np.where(a1 < -180., a1+360., a1)
        #a1 = np.where(a1 > +180., a1-360., a1)
        a1 = a1 * cos(0.5*(lat_v[j,i] + lat_v[j+1,i])*rad)
        angle_v = np.arctan2(a1, -a2)

        # Average
        angle0 = 0.5*(angle_u + angle_v)

        # pn, pm with difference formulation, using spherical distance
        pm0 = 1.0 / dist(lon_u[j,i+1], lat_u[j,i+1], lon_u[j,i], lat_u[j,i])
        pn0 = 1.0 / dist(lon_v[j+1,i], lat_v[j+1,i], lon_v[j,i], lat_v[j,i])

        # Central diffence formulation
        dndx0 = 0.5*(1/pn[j,i+1]-1/pn[j,i-1])
        dmde0 = 0.5*(1/pm[j+1,i]-1/pm[j-1,i])

        # Analytical expressons

        R = gmap.ellipsoid.a
        phi0 = gmap.lat_ts * rad
        m = gmap.map_scale(i,j)
        dx = gmap.dx
        x = (i - gmap.xp)*dx
        y = (j - gmap.yp)*dx
        A = - dx**2 / (m**2 * R**2 * (1+sin(phi0)))

        angle1 = (gmap.ylon - lon_rho[j,i])*rad
        pm1 = m / dx
        dndx1 = x * A
        dmde1 = y * A

        # Testing
        # Note: should work with variables saved as both
        #   float32 and float64 (only tested 64)
        self.assertAlmostEqual(angle[j,i], angle_u, places=2)
        self.assertAlmostEqual(angle[j,i], angle_v, places=2)
        self.assertAlmostEqual(angle[j,i], angle0, places=3)
        self.assertAlmostEqual(angle[j,i], angle1)

        self.assertAlmostEqual(pm[j,i], pm0, places=6)
        self.assertAlmostEqual(pn[j,i], pn0, places=6)
        self.assertAlmostEqual(pm[j,i], pm1)
        self.assertAlmostEqual(pn[j,i], pm1)

        self.assertAlmostEqual(dndx[j,i], dndx0)
        self.assertAlmostEqual(dmde[j,i], dmde0)
        self.assertAlmostEqual(dndx[j,i], dndx1, places=5)
        self.assertAlmostEqual(dmde[j,i], dmde1, places=5)

    # ---------------------------------------------------
         
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

