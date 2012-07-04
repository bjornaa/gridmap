=======================
Grid file specification
=======================

The ROMS grid file shall meet the following requirements

1. It should work without modifying ROMS

2. It should make clear what projection is being used

3. It should follow the CF-standard as far as possible

This is possible, since ROMS does not care about additional attributes
or variables, as long as it gets what it want.

[Nevne bruk av standard_names ....]


Global attributes
-----------------

ROMS require one global attribute, ``type``,  which can be anything.

The CF standard recommends ``history``, ``Conventions``, ``institution``.

In addition it can be useful to have a name for the grid, stored in a
global attribute ``gridname``.

The grid created by the ``initgrid.py`` example has::
 
  // global attributes:
       :gridname = "demo10km" ;
       :type = "ROMS grid file" ;
       :history = "Created by gridmap.define_grid.py" ;
       :Conventions = "CF-1.2" ;
       :institution = "Institute of Marine Research"


Dimensions
----------

As the grid is staggered, ROMS expect a several dimensions in the
:math:`\xi` and :math:`\eta` directions. They are named ``xi_m`` and
``eta_m`` where the modifier ``m`` indicates the stagerring, ``rho``,
``u``, ``v``, or ``psi``.

Not expected by ROMS, but common and useful in grid files is an extra
unlimited dimension ``bath``. This is used to index various more or less
raw bathimetries. This makes it possible to go back and smooth the
topography again, if needed.

Coordinate variables
--------------------

A `coordinate variable` is a 1D variable with the same name as a
dimension. The netcdf user's guilde and the CF-standard recommends the
use of such variables. Coordinate variables for the geographical
location of grid points are required as part of the grid mapping
definition. ROMS does not require coordinate variables.

As the polar stereographic grid is Cartesian in the projected plane,
coordinate variables like ``xi_rho`` is defined and gives the locations
in the projected plane in meters. Basically, the grid coordinates
multiplied by `dx`. They are defined as::
  
  float xi_rho(xi_rho) ;
      xi_rho:long_name = "X coordinate of RHO-points" ;
      xi_rho:standard_name = "projection_x_coordinate" ;
      xi_rho:units = "meter" ;

ROMS grid files sometimes has a series of 2D variables ``x_m`` and
``y_m``, with ``m`` in (``rho``, ``u``, ``v``, ``psi``), giving the locations of
the grid points in the projected plane. As these variables are not
required by ROMS or the CF-standard, they are omitted here. If needed
they are the Cartesian product of the coordinate variables. In numpy
notation::

  x_rho= np.meshgrid(eta_rho, xi_rho)  # Sjekk rekkefølge




Grid mapping variable
---------------------

A grid mapping variable is an empty variable containint the grid
mapping information as attributes. This is defined by the CF standard
and gives applications and users standardized information on the
projection used. The standard requires geographical coordinate
variables, longitude and latitude of the grid point

To follow the CF-standard new attributes are added to the data,
``standard_name`` when there is one, ``coordinates`` and
``mapping``::

  float h(eta_rho, xi_rho) ;
      h:long_name = "Final bathymetry at RHO-points" ;
      h:standard_name = "sea_floor_depth" ;
      h:units = "meter" ;
      h:field = "bath, scalar" ;
      h:coordinates = "lon_rho lat_rho" ;
      h:mapping = "grid_mapping" ;

The CF-standard http://cf-pcmdi.llnl.gov/ section 5.6 describes how to
indicate the horizontal coordinate system and the mapping projection. 
Version 1.2 of the standard gives conventions for the polar
stereographic mapping. This is done by defining a variable without
data and with the following attributes in the spherical case::

  int grid_mapping ;
      grid_mapping:long_name = "grid mapping" ;
      grid_mapping:grid_mapping_name = "polar_stereographic" ;
      grid_mapping:ellipsoid = "sphere" ;           // not required
      grid_mapping:earth_radius = 6371000. ;
      grid_mapping:latitude_of_projection_origin = 90. ; 
      grid_mapping:straight_vertical_longitude_from_pole = 58 ;  // ylon
      grid_mapping:standard_parallel = 60. ;                     // lat_ts
      grid_mapping:false_easting = 4182500. ;                    // xp*dx
      grid_mapping:false_northing = 2572500. ;                   // yp*dx 
      grid_mapping:dx = 10000 ;                     // not required
      grid_mapping:proj4string = "+proj=stere +R=6371000.0 +lat_0=90 +lat_ts=60.0 \\ 
                  +x_0=4182500.0 +y_0=2572500.0 +lon_0=58" ;

The name of the variable is not determined by the CF-standard. This
allows for several grids with different grid mappings in the same
NetCDF file. Here "grid_mapping" is suggested for polar stereographic
names. The name can be found by looking up the ``mapping`` attribute
of the mapped variable, for instance ``h`` as above.

The ``ellipsoid``, ``dx`` and ``proj4string`` attributes are not
required by the CF-standard, but considered useful. They are not
strictly necessary.  The grid spacing ``dx`` can found as the second
element in the ``xi_rho`` dimension variable and the proj4string can
be gathered from the other attributes.

For the ellipsoid, the default is a sphere. If a non-spherical
ellipsoid is used, it is identified by the value of the attributes
``semi_major_axis`` and ``inverse flattening`` replacing
``earth_radius``. For WGS84 this becomes::

  int grid_mapping ;
      grid_mapping:long_name = "grid mapping" ;
      grid_mapping:grid_mapping_name = "polar_stereographic" ;
      grid_mapping:ellipsoid = "WGS84" ;            // not required
      grid_mapping:semi_major_axis = 6378137. ;
      grid_mapping:inverse_flattening = 298.257223563 ;
      grid_mapping:latitude_of_projection_origin = 90. ;          
      grid_mapping:straight_vertical_longitude_from_pole = 58 ;  // ylon
      grid_mapping:standard_parallel = 60. ;                     // lat_ts
      grid_mapping:false_easting = 4182500. ;                    // xp*dx
      grid_mapping:false_northing = 2572500. ;                   // yp*dx
      grid_mapping:dx = 10000 ;                     // not required
      grid_mapping:proj4string = "+proj=stere +ellps=WGS84 +lat_0=90 \\
           +lat_ts=60.0 +x_0=4182500.0 +y_0=2572500.0 +lon_0=58" ;
  
Other variables
---------------

ROMS requires (but does not use) two variables ``xl`` and ``el``. They
give the domain size in meters. They can be anything, but slightly
useful values are ``Lm*dx`` and ``Mm*dx`` respectively. Also required
(and used) by ROMS is the character variable ``spherical``. This
should have the value "T" for *True*, even when used with a
non-spherical ellipsoid. The alternative is "F", which tells ROMS that
the grid is Cartesian and that the metrical terms should be neglected
[sjekk dette].
