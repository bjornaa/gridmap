=======================
Grid file specification
=======================

The ROMS grid file shall meet the following requirements

1. It should work without modifying ROMS

2. It should make clear what projection is being used

3. It should follow the CF-standard as far as possible

This is possible, since ROMS does not care about additional attributes
or variables, as long as it gets what it want.

Dimensions
----------

As the grid is staggered ROMS expect a series of dimensions in the
:math:`\xi` and :math:`\eta` directions. They are named `xi_m` and `eta_m`
where the modifier `m` indicates the stagerring, 'rho', 'u', 'v', or 'psi'.




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
