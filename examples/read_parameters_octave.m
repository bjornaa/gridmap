% read_parameters_octave.m
%
% Read polar stereographic parameters and grid size
% from a ROMS grid file (with polar stereographic info)
% using octave with the octcdf toolbox
% 
% ---------------------------------
% Bjørn Ådlandsvik <bjorn@imr.no>
% Institute of Marine Research
% 2012-07-17
% ---------------------------------

gridfile = 'demo10km_grd.nc';

nc = netcdf(gridfile, 'r');

Lm = length(nc('xi_rho')) - 2;
Mm = length(nc('eta_rho')) - 2;

v = nc{'grid_mapping'};
dx = v.dx;
xp = v.false_easting / dx;
yp = v.false_northing / dx;
ylon = v.straight_vertical_longitude_from_pole;

printf("%s\n", gridfile)
printf("Lm, Mm = %d, %d\n", Lm, Mm)
printf("xp, yp, dx, ylon = %f, %f, %f, %f\n", xp, yp, dx, ylon)

