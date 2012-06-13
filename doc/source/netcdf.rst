=======================
Grid file specification
=======================

The grid file follows the CF-standard http://cf-pcmdi.llnl.gov/.
Version 1.2 of the standard gives conventions for the polar
stereographic mapping. This is done by defining a variable without
data and attributes::

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
  false_easting = `xp`
  false_northing = yp

In addition the following general mapping attributes are used
in the sphere case::

  earth_radius = 6371000

and for WGS84::

  semi_major_axis = 6378137.0
  inverse_flattening = 298.257223563
  
Two extra attributes not required by the CF-standard are
useful additions::

  dx = `dx`
  proj4string = option string used by proj4 to recreate the mapping

Spørsmål: Er det mer korrekt å bruke.

  false_easting = dx * xp
  
Gjør det i proj4string. Da blir alt i meter.
