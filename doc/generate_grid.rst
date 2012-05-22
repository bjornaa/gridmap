==========================================
Howto: Generate a polar stereographic grid
==========================================

:Author: Bjørn Ådlandsvik
:Organization: Institute of Marine Research
:Contact: bjorn@imr.no
:Date: 2012-04-25

This document shows a procedure for generating a polar stereographic
grid usin python tools. The method is modular, so that the separate
tasks can be done by other tools if they are better suited.

The document is written in restructured text and can easily be
converted to html or pdf (via latex).

Define the grid and initialize the grid file
--------------------------------------------

This step defines the grid and initializes the grid file.

The grid is given by parameters::

  Lm = The number of interior grid cells in X-direction
  Mm = The number of interior grid cells in Y-direction
 
  xp = X grid coordinate of the North Pole
  yp = Y grid coordinate of the North Pole
  dx = grid resolution (at latitude of true scale)
  ylon = vertical meridian
  
  lon_ts = longitude of true scale, always 60.0
  ellipsoid: sphere (with radius) or WGS84
   
Let X and Y denote grid coordinates i.e.
X = i, Y = j is the (i,j)-th rho-point with ROMS indexing
(starting from zero for rho-points).

Adjust the grid definition variables in the user setting
section of the script define_grid.py

The script generates a grid file with place-holders for all variables
and fills in topography-independent values, such as longitudes,
latitudes, metric factors, coriolis and rotation angle.

The script adds metainformation including a map projection details
in a format complaying with the CF-standard (CF-1.2).

The grid with coastline can be plotted by plot_basegrid.py (requiring
basemap) or plot_base.. (requirs GMT).

Add the raw topography
----------------------

If the topography comes from a finer source than the target grid
the method of the ``etopo1grid.py`` can be used. As there are many 
topography values inside each grid cell, it returns the average of the
sea values inside the grid cell. If the grid cell is entirely at land
a value of zero is used. This is added to the grid file as the first
record in the hraw variable. 

The script can be speeded up by restricting the range of longitude
and latitude to be read in. The script may take a long time, in
particular if the North Pole is included in the grid.
(45 minutes for the arctic10_WGS84 grid). It is also recommended to
make a backup copy of the grid file at this stage, in case further
processing ruins some of the information.

The script also makes an auxiliary netcdf file, ``aux.nc``, saving the
mean value and standard deviation of topography within each grid cell
as well as the total number C0 of topography points and the number of
sea points C inside the cell.

If the resolution of the topography source is similar to the target
grid, the method above is not feasible. Instead an interpolation
routine must be used.

Extract a coast line
--------------------

The script ``makecoast.py`` uses basemap to make a coast line. In the user
settings the name of the grid file and the GSHHS resolution code must
be given. An area threshold (default to half a grid cell area) is also
used to avoid many small islands. 

The resulting coast line is given in grid coordinates and is clipped
towards the grid boundary (the outer boundary of the boundary cells).

NOTE: A good idea to add the coast line to the grid file?



