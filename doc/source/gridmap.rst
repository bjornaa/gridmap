==============================
The gridmap package for python
==============================

The gridmap package defines a class `PolarStereographic`.
A typical usage is::

  import gridmap
  ...
  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)

this sets up the projection. If the grid size is needed
(for plotting or creating a grid file) the number of internal grid
cells, ``Lm``, ``Mm`` are optional arguments. Other optional arguments 
are the ellipsoid (defalt a sphere of radius 6371 km) and the 
latitude of true scale (default is 60°N) 
The full interface is::

  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, Lm=None, Mm=None,
                      ellipsoid=gridmap.sphere, lat_ts=60.0)

:mod:`gridmap`

.. automodule:: gridmap
   :members:
