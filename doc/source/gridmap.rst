==============================
The gridmap package for python
==============================

The gridmap package defines a class `PolarStereographic`.
A typical usage is::

  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)

this sets up the projection. If the grid size is needed
(say for plotting) it can be added as a two-tuple (Lm, Mm).
The WGS84 ellipsoid is available by an optional argument::

  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, 
       	 			    ellipsoid=gridmap.WGS84)


