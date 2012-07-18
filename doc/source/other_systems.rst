=======================================================
Handling Polar Stereographic grids in different systems
=======================================================

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
  dx = v.dx
  xp = v.false_easting / dx
  yp = v.false_northing / dx
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

Using the fortran 90 

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
stereographic parametere from the grid mapping variable

  ncdump -h $GRIDFILE > $CDL_FILE
  ...
  line=( `grep "xi_rho = " $CDL_FILE` )
  Lp=${line[2]}
  let Lm=Lp-2
  ...
  line=( `grep grid_mapping:straight_vertical_longitude_from_pole $CDL_FILE` )
  ylon=${line[2]}

Octave
------

Octave (http://www.gnu.org/software/octave/) is an open source
Matlab-clone. There is a NetCDF interface, ``octcdf`` available
from http://octave.sourceforge.net/octcdf/index.html .

::

  nc = netcdf(gridfile, 'r');
  ...
  Lm = length(nc('xi_rho')) - 2;
  ...
  v = nc{'grid_mapping'};
  dx = v.dx;
  xp = v.false_easting / dx;
  yp = v.false_northing / dx;
  ylon = v.straight_vertical_longitude_from_pole;






Matlab
------

Newer versions of Matlab comes with a NetCDF interface. 
[Har ikke Matlab-lisens, noen bør teste] 





Computing grid coordinates
==========================

Here a standardized example is used. A grid mapping for the Arctic is
defined by ``xp`` = 418.25, ``yp`` = 257.25, ``dx`` = 10000, and
``ylon`` = 58 on a spherical earth (radius=6371km). These parameters can
be read from the file ``demo10km_grid.nc`` as above or given
explicitly. The projection does not depend on the grid size, so
``Lm``, ``Mm`` are left undefined.

+--------+--------+-------+------+
| xp     |    yp  |   dx  | ylon |
+========+========+=======+======+
| 418.25 | 257.25 | 10000 | 58   |
+--------+--------+-------+------+

Litt tekst


xp       yp       dx     ylon
======   ======   =====  ====
418.25   257.25   10000  58
======   ======   =====  ====

Litt tekst

  +--------+-----+-----+---------------+---------------+
  |        | lon | lat | x             |  y            |
  +========+=====+=====+===============+===============+
  | sphere | 2   | 66  | 208.754891658 | 115.943765186 |
  +--------+-----+-----+---------------+---------------+
  | WGS84  | 2   | 66  | 207.924459414 | 115.383631565 |
  +--------+-----+-----+---------------+---------------+



The different systems are used to compute the grid coordinates of
station M (2°E, 66°N). The correct values (as computed by proj4) are::

  x = 208.754891658 
  y = 115.943765186 

The `inverse` or `backwards` projection is used to find the location
of the origin (rho-point in lower left grid cell, x=0, y=0). Here the
authorative values from proj4 are

   +--------+-----+-----+-----------------+---------------+
   |        | x   | y   | lon             | lat           |
   +========+=====+=====+=================+===============+
   | sphere | 0.0 | 0.0 | -0.405875241137 | 45.1156804622 |
   +--------+-----+-----+-----------------+---------------+
   | WGS84  | 0.0 | 0.0 | -0.405875241137 | 45.2201521896 |
   +--------+-----+-----+-----------------+---------------+


::

  lon = -0.405875241137
  lat = 45.2201521896


python with gridmap
-------------------

The ``Gridmap`` package for python is partly designed for this
purpose. The code is::

  import gridmap
  ...
  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)
  ...
  lon, lat = 2, 66
  x, y = gmap.ll2grid(lon, lat)
  ...
  x, y = 0, 0
  lon, lat = gmap.grid2ll(x, y)

With the numbers of decimals above, the results are identical to the
proj4 control values.

proj4
-----

Proj4 http://trac.osgeo.org/proj/ is perhaps the standard software for
map projections offering a wealth of projections and ellipsoids.

proj4 from shell
................

Proj4 comes with a command line program ``proj`` that can be used from
the shell to do projections. A polar stereographic mapping with values
for ``XP``, ``YP``, ``DX``, ``YLON`` on a spherical earth is specified
by::

  XPDX=$(echo "($XP*$DX)" | bc)   # multiply floats XP and DX
  YPDX=$(echo "($YP*$DX)" | bc)
  proj -m 1:$DX +proj=stere +R=6371000 +lat_0=90 +lat_ts=60 +x_0=$XPDX +y_0=$YPDX +lon_0=$YLON

Here the option ``-m 1:$DX`` scales the output to directly to grid
coordinates. In the WGS84 case, the option ``+R=637100`` is replaced by 
``+ellps=WGS84``. 

The input is taken from a file, or as done here, a ``here`` construct
in the shell.

The inverse projection is done by the option ``-I`` or the command 
``invproj`` with the same options. 

When a grid file is available, the ``proj4string`` simplifies the
script, providing the options and eliminates the complicated
multiplication of floating points numbers in the shell. The usage is
then simply::

  PROJ4STRING=`ncdump -h $GRIDFILE | grep proj4string`
  $PROJ $PROJ4STRING

Note, as the option ``-m`` is not used, the results are
unscled and has to be divided by ``DX``.
[test at dette virker]



proj4 from python
-----------------

Proj4 can be run from a python script instead of the shell. This
offers some advantages. It is system independent and can even be used
on systems like Microsoft Windows which does not come with a standard
shell. [Test dette] It uses a real programming language and can be
integrated with other packages like ``gridmap``. The ``subprocess``
module from the python standard library gives complete control of the
``proj`` program.

The example script ``project_proj4.py`` provides two functions
``proj`` and ``invproj`` handling the subprocess details. The
``PolarStereographic`` class in ``gridmap`` has an attribute
``proj4string`` which contains the projection options (the ones with
plusses). The use is as simple as::

  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)
  x, y = proj(gmap.proj4string, lon, lat)
  x, y = x/dx, y/dx

Note that the ``proj`` function is meant to be general for the use of proj4,
not only for polar stereographic grids it does not do the scaling to
grid coordinates. 


GMT
---

The generic mapping tools ``GMT`` (ref) is a much used package for map
projections and plotting. It is independent of proj4 and implements
the map projections separately. It is usually used from the shell, but
similarly to proj4, usage from python is recommended. Also for GMT the
a complicated sequence of command line options is required.

The example script ``project_gmt.py`` introduces projection functions
``proj`` and ``invproj`` [faktisk kalt ``mapproject`` using the GMT program ``mapproject`` similar
to the proj4-functions. An example usage::

  [Virker ikke for spherical earth???, feil i ellipsoid???]


  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)
  projection = '-Js%s/90.0/%s/1:%s' % \
                 (str(ylon), str(gmap.lat_ts), str(100*dx))
  ellipsoid = "--ELLIPSOID=%s" % str(gmap.ellipsoid.a)
  extent = '-R0/1/60/61'   # Map extent, actual values are not used
  offset = '-C%s/%s' % (str(xp), str(yp))
  gmtstring = " ".join((ellipsoid, projection, offset, extent))

  x, y = mapproject(gmtstring, lon, lat)

This gives grid coordinates directly, as the scaling is part of the
``-J``-string.

[Lage en GMTstring i gridmap-pakke??]

Basemap
-------

``Basemap`` is a mapping toolkit for the python package matplotlib. It
is developed by Jeffrey Whitaker with web site
http://matplotlib.github.com/basemap/. The map projection part is
based on ``proj4``. The projection can be done by putting the origin
at the North Pole (with arbitraty values for ``urcrnrlon`` and
``urcrnrlat``). After rescaling ``xp`` and ``yp`` gives the required
offset::

  from mpl_toolkits.basemap import Basemap
  ..

  # Set up the projectiob
  bmap = Basemap(projection='stere', rsphere=6371000,
                 llcrnrlon = ylon, llcrnrlat = 90.0,
                 urcrnrlon = ylon+1, urcrnrlat = 89.0,  
                 lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)
  
  # Do the projection
  x, y = bmap(lon, lat) 
  x = x / dx + xp   # rescale and add offset to get grid coordinates
  y = y / dx + yo

For ``WGS84`` use::

  # Elliptical parameters for WGS84
  a = 6378137.0
  f = 1./298.257223563
  b = a*(1-f) 

  # Projection in basemap 
  bmap = Basemap(projection='stere', rsphere=(a, b),
            llcrnrlon = ylon, llcrnrlat = 90.0,
            urcrnrlon = ylon+1, urcrnrlat = 89.0,  
            lat_0 = 90.0, lon_0 = ylon, lat_ts = 60.0)

The ``gridmap`` package has support for basemap, simplifying the
``Basemap`` interface. The origin is moved to the grid origin,
eliminating the need for offset::

  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)
  bmap = gmap.basemap()
  x, y = bmap(lon, lat) 
  x = x / dx # rescale to get grid coordinates
  y = y / dx 

  # OOPS, trenger gmap.Lm, gmap.Mm, fiks dette


pyproj
------

Jeffrey Whitaker (staving korrekt) has implementet a interface to
``proj4`` as library.  It is part of ``basemap``, but also available
separately https://code.google.com/p/pyproj/. This is the recommended
way to use ``proj4`` from ``python``, as it does not need the
``subprocess`` package and can be used effectively on ``numpy``
arrays.

As a separate package it is imported as::

  import pyproj

From ``basemap`` it is imported as::

  import mpl_toolkits.basemap.pyproj as pyproj
  (sjekk dette)

The user interface is simplest by using ``gridmap`` and 
the ``proj4string``::

  gmap = gridmap.PolarStereographic(xp, yp, dx, ylon)
  pmap = pyproj.Proj(gmap.proj4string)
  x, y = pmap(lon, lat)
  x, y = x/dx, y/dx

  # The inverse projection
  lon, lat = pmap(x*dx, y*dx, inverse=True)


M-Map
-----

M_Map is a maping package for Matlab and Octave. 


FIMEX
-----



map plots
=========





