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

The workflow is modular and incremental, leaving it open for different
methods.  After the grid file is initiated in step 2, other tools can
be used to add topography, edit the landmask and smooth the
topography. This can be tools in fortran, matlab or python.


1. Define the grid
------------------

A polar stereographic grid is given by the following parameters::

  Lm = The number of interior grid cells in X-direction
  Mm = The number of interior grid cells in Y-direction
 
  xp = X grid coordinate of the North Pole
  yp = Y grid coordinate of the North Pole
  dx = grid resolution (at latitude of true scale)
  ylon = vertical meridian
  
  lon_ts = longitude of true scale, always 60.0
  ellipsoid: sphere with radius 6371 km (default) or WGS84
   
The first step is to determine suitable values here. Some tools are
available.

The least intuitive is to find `xp` and `yp`. If you have decided 
the resolution `dx` and rotation angle `ylon` and know where to put
the origin. The script `examples/find_xpyp.py` may be useful. Copy the 
script to the working directory, edit the 'user settings' and run::

  python find_xpyp

Choose suitably round values for `xp` and `yp` close to the values returned by
the script.

The next step is to plot the domain. If you have the `basemap` toolkit
for `matplotlib` installed, copy `examples/plot_basegrid.py` to the
working directory, edit the parameters in the file and run::

  python plot_basegrid.py

For small adjustements, increasing `xp` resp. `yp` extends the domain
to the left resp. downwards. Increasing `Lm` and `Mm` extends the
domain to the right resp. upwards.

If the Generic Mapping Tools (GMT) are installed on the system, the script
`examples/plot_gmt_map.py` can be used similarly.

[legg inn matlab/python plot med mmap]

2. Initiate the grid file
-------------------------

Copy the script `examples/initgrid.py` to your working directory and 
edit the mapping parameters, grid name and institution.  Run::

  python initgrid.py

The script generates a grid file with place-holders for all variables
and fills in topography-independent values, such as longitudes,
latitudes, metric factors, Coriolis and rotation angle.

The script adds metainformation including a map projection details
in a format complying with the CF-standard (CF-1.2).

3. Add the raw topography
-------------------------

This step is quite open, as there are a lot of different bathymetric
datasets available. If the resolution in the data source is similar to
the grid use your favourite interpolation method. If the data source
is finer, the averaging procedure described here is recommended.
If the data source is coarser than the grid, you have a problem.

In any case, save the unmodified topography as the first field in
`hraw` to make it possible to go back and redo the later steps if
wanted. It is also recommended to make a backup copy of the file
after this step, as further modifications of the file may ruin it.

In the lucky case where a fine dataset is available, averaging the
topography points inside a grid cell is recommended. First, this is
simpler than interpolation and second it gives a better estimate of
the mean depth of the cell. The method is used in the script
`etopo1grid.py`. The script may take a long time to run, in particular
if the North Pole is included in the grid. Restricting the range of 
topography points considered can speed it up substantially.

The script also makes an auxiliary netcdf file, ``aux.nc``, saving the
mean value and standard deviation of topography within each grid cell
as well as the total number C0 of topography points and the number of
sea points C inside the cell. This file may be useful for further
steps in the grid generation.

4. Extract a coast line
-----------------------

The script ``makecoast.py`` uses basemap to make a coast line. In the user
settings the name of the grid file and the GSHHS resolution code must
be given. An area threshold (default to half a grid cell area) is also
used to avoid many small islands. 

The resulting coast line is given in grid coordinates and is clipped
towards the grid boundary (the outer boundary of the boundary cells).

NOTE: A good idea to add the coast line to the grid file?

If the matlab tools are used, this is not necessary as a separate
tool, as `editmask.m` makes a coast line in grid coordinates the first
time it is called.

5. Make an initial land mask
----------------------------

6. Edit the land mask
---------------------

7. Compute the other masks
--------------------------

8. Smooth the topography
------------------------

This step may be quite demanding. The topography should be smoothed to
1) set a minimum depth, 2) avoid steep gradients causing problems with
the pressure gradient and 3) avoid noice on the smallest (2*dx) wave
length.





