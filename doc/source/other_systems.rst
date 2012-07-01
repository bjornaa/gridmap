========================================
Polarsterographic grids in other systems
========================================

Reading the parameters from file
================================

To be consistent and avoid errors scripts should not hardcode the
polar stereographic parameters. Instead this information should be
read from the grid file. The information is coded as attributes of a
variable. As a netCDF file may contain variables with different grid
mappings, this name of the variable is not specified by the
CF-standard. It can be found from the attribute
``grid_mapping`` of the ``h``-variable. To simplify, it is assumed
here that the variable is simply called ``grid_mapping``

Python with netCDF4-python
--------------------------

Without the ``gridmap`` package this can be done in a straight-forward way:: 

  from netCDF4 import Dataset
  ...
  f = Dataset(grid_file)
  v = f.variables['grid_mapping']
  xp = v.false_easting / v.dx
  yp = v.false_northing / v.dx
  dx = v.dx
  ylon = v.straight_vertical_longitude_from_pole
  Mm = len(f.dimensions['xi_rho'])-2
  Lm = len(f.dimensions['eta_rho'])-2

With the ``gridmap`` package, it is simplified to::

  from netCDF4 import Dataset
  import gridmap
  ...
  f = Dataset(grid_file)
  gmap = gridmap.fromfile(grid_file)

Fortran
-------

::

  use netcdf
  ...
  integer :: status, ncid, dimid, varid
  ...
  status = nf90_open(gridfile, nf90_nowrite, ncid)
  ...
  status = nf90_inq_dimid(ncid, 'xi_rho', dimid)
  status = nf90_inquire_dimension(ncid, dimid, len=Lm)
  Lm = Lm - 2
  ...
  status = nf90_inq_varid(ncid, 'grid_mapping', varid)
  status = nf90_get_att(ncid, varid, 'dx', dx)



shell
-----

The ``bash``-shell does not have a shell NetCDF reader. But the
information can be gained from the CDL-file created by ``ncdump``.
The grid size can be found from the dimension settings and the
stereographic parametere from the grid mapping variable::

  ncdump -h $GRIDFILE > $CDL_FILE
  ...
  line=( `grep "xi_rho = " $CDL_FILE` )
  Lp=${line[2]}
  let Lm=Lp-2
  ...
  line=( `grep grid_mapping:straight_vertical_longitude_from_pole $CDL_FILE` )
  ylon=${line[2]}







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


