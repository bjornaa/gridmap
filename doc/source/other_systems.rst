========================================
Polarsterographic grids in other systems
========================================


proj4
=====

Proj4 http://trac.osgeo.org/proj/ is perhaps the standard software for
map projections offering a wealth of projections and ellipsoids.

A mapping with values for XP, YP, DX, YLON can be 
be recreated in a shell script by::

  XPDX=$(echo "($XP*$DX)" | bc)
  YPDX=$(echo "($YP*$DX)" | bc)
  proj -m 1:$DX +proj=stere +R=6371000 +lat_0=90 +lat_ts=60 +x_0=$XPDX +y_0=$YPDX +lon_0=$YLON

This is for a spherical earth. For the WGS84 ellipsoid use::

  proj -m 1:$DX +proj=stere +ellps=WGS84 +lat_0=90 +lat_ts=60 +x_0=$XPDX +y_0=$YPDX +lon_0=$YLON

An example can be found in the shell script `project_proj4.sh`

proj4string
-----------

The use of proj is greatly simplified by projstrings containing all
the arguments. Such a string is
provided in the ROMS grid files adhering the standard here
and can be found by::
 
  ncdump -h $GRIDFILE | grep projstring

Using python with the gridmap package, the projstring can be found
by::

  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon, ...)
  gmap.projstring

GMT
===

Basemap
=======

Basemap is a mapping package for the python package matplotlib. It is
developed by Jeffrey Whitaker with web site
http://matplotlib.github.com/basemap/. The map projection part is
based on proj4.

projections
-----------

map plots
---------

pyproj
======

Is a python interface to the proj4 library. It is part of basemap,
but also available separately https://code.google.com/p/pyproj/

M-Map
=====

M_Map is a maping package for Matlab and Octave. 


FIMEX
=====


