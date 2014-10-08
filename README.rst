GRIDMAP
=======

Python toolbox for ROMS grid files based on explicit map projections.

Many approaches, such as `Seagrid` (Matlab), `Gridpak` (Fortran),
`Gridgen` (python) are using advanced numerics to generate orhogonal
grids based on prescribed boundaries. This gives great flexibility for
adapting the grid to the coastal geometry. The downside is that there
is no obvious mapping between longitude/latitude and the grid
coordinates.

Gridmap has a different approach. The grid is based on well documented
orhogonal map projections. The generated grid is therefore
automatically orhogonal (in most cases even orthonormal) and there is
an explicit transformation between lon/lat and grid coordinates.  It
is presently independent of solvers in Fortran or C making it easy to
install. At some point it may develop to depend on the `pyproj`
package to benefit from the large range of projections provided.
Presently it focuses on the polar stereographic projection, a good
choice at relative high latitudes.

Bjørn Ådlandsvik  <bjorn at imr.no>
Institute of Marine Research
2014-10-08
