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

ROMS require one global attribute, `type`,  which can be anything.
A reasonable alternative "ROMS grid file".

The CF standard recommends `history`, `Conventions`, `institution`.

In addition the attribute `gridname` can be useful.

The grid created by the initgrid.py example has::
 
  // global attributes:
       :gridname = "demo10km" ;
       :type = "ROMS grid file" ;
       :history = "Created by gridmap.define_grid.py" ;
       :Conventions = "CF-1.2" ;
       :institution = "Institute of Marine Research"


Dimensions
----------

As the grid is staggered ROMS expect a several dimensions in the
:math:`\xi` and :math:`\eta` directions. They are named `xi_m` and `eta_m`
where the modifier `m` indicates the stagerring, 'rho', 'u', 'v', or
'psi'.

Not expected by ROMS, but common and useful in grid files is an extra
unlimited dimension `bath`. This is used to index various more or less
raw bathimetries.

Coordinate variables
--------------------

A `coordinate variable` is a 1D variable with the same name as a
dimension. The netcdf user's guilde and the CF-standard recommends the
use of such variables. Coordinate variables for the geographical
location of grid points are required as part of the grid mapping
definition. ROMS does not require coordinate variables.

As the polar stereographic grid is Cartesian in the projected plane,
coordinate variables like `xi_rho` is defined and gives the locations
in the projected plane in meters. Basically, the grid coordinates
multiplied by `dx`. They are defined as::
  
  float xi_rho(xi_rho) ;
       xi_rho:long_name = "X coordinate of RHO-points" ;
       xi_rho:standard_name = "projection_x_coordinate" ;
       xi_rho:units = "meter" ;

ROMS grid files sometimes has a series of 2D variables `x_m` and
`y_m`, with 'm' in ('rho', 'u', 'v', 'psi'), giving the locations of
the grid points in the projected plane. As these variables are not
required by ROMS or the CF-standard, they are omitted here. If needed
they are the Cartesian product of the coordinate variables. In numpy
notation::

  x_rho= np.meshgrid(eta_rho, xi_rho)  # Sjekk rekkefølge

ROMS require (but do not use) two variables `xl` and `el`. They give
the domain size in meters. They can be anything, but slightly useful
values are `Lm*dx` and `Mm*dx` respectively.


Grid mapping variable
---------------------

A grid mapping variable is an empty variable containint the
grid mapping information as attributes. This is defined by the CF
standard and gives applications and users standardized information
on the projection used.

For a coor

The standard requires geographical coordinate variables, longitude and
latitude of the grid point

To follow the CF-standard new attributes are added to the data,
`standard_name` when there is one, `coordinates` and
`mapping`::

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
data and with the following attributes::

  grid_mapping_name
  straight_vertical_longitude_from_pole
  latitude_of_projection_origin - Either +90. or -90.
  Either standard_parallel or scale_factor_at_projection_origin
  false_easting
  false_northing

If the grid is defined by `xp`, `yp`, `dx`, `ylon` as usual, these 
parametes are set to:: 

  grid_mapping_name = "polar_stereographic"
  straight_vertical_longitude_from_pole = ylon
  1atitude_of_projection origin = 90.0
  standard_parallell = 60.0
  false_easting = xpx
  false_northing = yp

In addition the following general mapping attributes are used
in the sphere case::

  earth_radius = 6371000

and for WGS84::

  semi_major_axis = 6378137.0
  inverse_flattening = 298.257223563
  
Two extra attributes not required by the CF-standard are
useful additions::

  dx = dx
  proj4string = option string used by proj4 to recreate the mapping

Spørsmål: Er det mer korrekt å bruke.

  false_easting = dx * xp
  
Gjør det i proj4string. Da blir alt i meter.
