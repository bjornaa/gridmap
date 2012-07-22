program read_parameters

! Fortran 90 program for reading polar stereographic parameters
! from a ROMS grid file (with this information)

! On my laptop it compiles with:
! gfortran -I/usr/include read_parameters.f90 -l netcdff

! ----------------------------------------
! Bjørn Ådlandsvik <bjorn@imr.no>
! Institute of Marine Research
! 2012-06-30
! ----------------------------------------

use netcdf
 
character(len=*), parameter :: gridfile = "demo10km_grid.nc"

real :: xp, yp, dx, ylon
integer :: Lm, Mm

integer :: status, ncid, dimid, varid


status = nf90_open(gridfile, nf90_nowrite, ncid)

status = nf90_inq_dimid(ncid, 'xi_rho', dimid)
status = nf90_inquire_dimension(ncid, dimid, len=Lm)
Lm = Lm - 2

status = nf90_inq_dimid(ncid, 'eta_rho', dimid)
status = nf90_inquire_dimension(ncid, dimid, len=Mm)
Mm = Mm - 2

status = nf90_inq_varid(ncid, 'grid_mapping', varid)
status = nf90_get_att(ncid, varid, 'dx', dx)

status = nf90_inq_varid(ncid, 'grid_mapping', varid)
status = nf90_get_att(ncid, varid, 'false_easting', xp)
xp = xp / dx

status = nf90_inq_varid(ncid, 'grid_mapping', varid)
status = nf90_get_att(ncid, varid, 'false_northing', yp)
yp = yp / dx

status = nf90_inq_varid(ncid, 'grid_mapping', varid)
status = nf90_get_att(ncid, varid,                           &
           'straight_vertical_longitude_from_pole', ylon)


print *, "Lm, Mm =", Lm, Mm
print *, "xp, yp, dx, ylon = ", xp, yp, dx, ylon



end program read_parameters

