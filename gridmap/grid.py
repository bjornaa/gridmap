# -*- coding: utf-8 -*-

"""Define a ROMS grid file, and fill in the map terms

CF compliant version for polarstereographic grid
on sphere or WGS84 ellipsoid

"""

# -------------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-04-15
# -------------------------------------

import numpy as np
from netCDF4 import Dataset
#from gridmap import fromfile
import gridmap  # fiks relativ import


# TODO: Make it possible to specify the order of the
#       global attributes
def make_empty_gridfile(grid_name, file_name, Lm, Mm, 
                        global_attributes, format):
    """
    Create a new empty ROMS grid file

    Arguments:
      gmap      : a gridmap.PolarStereographic instance
      grid_name : name of the grid
      file_name : name of the grid file,
                  default = '' giving grid_name + '_grid.nc'
      format    : 'NETCDF3_CLASSIC' or 'NETCDF4_CLASSIC'
                  default = 'NETCDF3_CLASSIC'

    Make space for all variables, including coordinate variables
    
    """
    
    if not file_name:  # Use default
        file_name = grid_name + '_grid.nc'

    # Ta denne med polarstereografisk ???
    gridmap_varname = 'grid_mapping' # Name of grid mapping variable
    
    # -----------------------
    # NetCDF file definition
    # -----------------------

    nc = Dataset(file_name, 'w', format=format)

    # Global attributes
    # Defaults
    global_defaults = dict(gridname    = grid_name, 
                           type        = 'ROMS grid file',
                           history     = 'Created by gridmap', 
                           Conventions = 'CF-1.2')
    d = {}
    d.update(global_defaults, **global_attributes)

    #print d
    for att, value in d.iteritems():
        setattr(nc, att, value)

    # Dimensions

    L,  M  = Lm+1, Mm+1
    Lp, Mp = Lm+2, Mm+2
    nc.createDimension('xi_rho',  Lp)
    nc.createDimension('eta_rho', Mp)
    nc.createDimension('xi_u',    L)
    nc.createDimension('eta_u',   Mp)
    nc.createDimension('xi_v',    Lp)
    nc.createDimension('eta_v',   M)
    nc.createDimension('xi_psi',  L)
    nc.createDimension('eta_psi', M)
    nc.createDimension('bath',    None)

    # --- Coordinate variables --- 
    # Not required by ROMS, but recommended by the CF standard

    v = nc.createVariable('xi_rho', 'd', ('xi_rho',), zlib=True)
    v.long_name = "X coordinate of RHO-points"
    v.standard_name = "projection_x_coordinate"
    v.units = "meter"

    v = nc.createVariable('eta_rho', 'd', ('eta_rho',), zlib=True)
    v.long_name = "Y coordinate of RHO-points"
    v.standard_name = "projection_y_coordinate"
    v.units = "meter"

    v = nc.createVariable('xi_u', 'd', ('xi_u',), zlib=True)
    v.long_name = "X coordinate of U-points"
    v.standard_name = "projection_x_coordinate"
    v.units = "meter"

    v = nc.createVariable('eta_u', 'd', ('eta_u',), zlib=True)
    v.long_name = "Y coordinate of U-points"
    v.standard_name = "projection_y_coordinate"
    v.units = "meter"

    v = nc.createVariable('xi_v', 'd', ('xi_v',), zlib=True)
    v.long_name = "X coordinate of V-points"
    v.standard_name = "projection_x_coordinate"
    v.units = "meter"

    v = nc.createVariable('eta_v', 'd', ('eta_v',), zlib=True)
    v.long_name = "Y coordinate of V-points"
    v.standard_name = "projection_y_coordinate"
    v.units = "meter"

    # --- Geographic variables

    v = nc.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "longitude of RHO-points"
    v.standard_name = "longitude"
    v.units = "degrees_east"

    v = nc.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "latitude of RHO-points"
    v.standard_name = "latitude"
    v.units = "degrees_north"

    v = nc.createVariable('lon_u', 'd', ('eta_u', 'xi_u'), zlib=True)
    v.long_name = "longitude of U-points"
    v.standard_name = "longitude"
    v.units = "degrees_east"

    v = nc.createVariable('lat_u', 'd', ('eta_u', 'xi_u'), zlib=True)
    v.long_name = "latitude of U-points"
    v.standard_name = "latitude"
    v.units = "degrees_north"

    v = nc.createVariable('lon_v', 'd', ('eta_v', 'xi_v'), zlib=True)
    v.long_name = "longitude of V-points"
    v.standard_name = "longitude"
    v.units = "degrees_east"

    v = nc.createVariable('lat_v', 'd', ('eta_v', 'xi_v'), zlib=True)
    v.long_name = "latitude of V-points"
    v.standard_name = "latitude"
    v.units = "degrees_north"

    v = nc.createVariable('angle', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "angle between xi axis and east"
    v.standard_name = "angle_of_rotation_from_east_to_x"
    v.units = "radian"
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    # Note: use spherical even if WGS84 ellipsoid
    # the alternative is cartesian without metric terms
    v = nc.createVariable('spherical', 'c', ())
    v.long_name = "grid type logical switch" 
    v.flag_values = "T, F"
    v.flag_meanings = "spherical Cartesian"

    # --- Metric terms

    v = nc.createVariable('pm', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "curvilinear coordinate metric in XI"
    v.units = "meter-1"
    v.field = "pm, scalar"
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    v = nc.createVariable('pn', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "curvilinear coordinate metric in ETA"
    v.units = "meter-1"
    v.field = "pn, scalar"
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    v = nc.createVariable('dndx', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "xi derivative of inverse metric factor pn"
    v.units = "meter"
    v.field = "dndx, scalar"
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    v = nc.createVariable('dmde', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "eta derivative of inverse metric factor pm"
    v.units = "meter"
    v.field = "pn, scalar"
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    # --- Grid map

    v = nc.createVariable(gridmap_varname, 'i', ())
    v.long_name = "grid mapping"
    #d = gmap.CFmapping_dict()
    #for att in d:
    #    setattr(v, att, d[att])
    #v.proj4string = gmap.proj4string

    # --- Topography

    v = nc.createVariable('hraw', 'd', ('bath', 'eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "Working bathymetry at RHO-points"
    v.standard_name = "sea_floor_depth"
    v.units = "meter" 
    v.field = "bath, scalar"
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    v = nc.createVariable('h', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "Final bathymetry at RHO-points"
    v.standard_name = "sea_floor_depth"
    v.units = "meter"
    v.field = "bath, scalar"
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    # --- Coriolis

    v = nc.createVariable('f', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = 'Coriolis parameter at RHO-points'
    v.standard_name = "coriolis_parameter"
    v.units = 'second-1'
    v.field = 'Coriolis, scalar'
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    # --- Masks 

    v = nc.createVariable('mask_rho', 'd', ('eta_rho', 'xi_rho'), zlib=True)
    v.long_name = "mask on RHO-points" 
    #v.standard_name = "sea_binary_mask"   # Not in standard table 
    v.option_0 = "land" 
    v.option_1 = "water" 
    v.coordinates = "lon_rho lat_rho"
    v.grid_mapping = gridmap_varname

    v = nc.createVariable ('mask_u', 'd', ('eta_u', 'xi_u'), zlib=True)
    v.long_name = "mask on U-points" 
    v.option_0 = "land" 
    v.option_1 = "water" 
    v.coordinates = "lon_u lat_u"
    v.grid_mapping = gridmap_varname

    v = nc.createVariable('mask_v', 'd', ('eta_v', 'xi_v'), zlib=True)
    v.long_name = "mask on V-points" 
    v.option_0 = "land" 
    v.option_1 = "water" 
    v.coordinates = "lon_v lat_v"
    v.grid_mapping = gridmap_varname

    v = nc.createVariable('mask_psi', 'd', ('eta_psi', 'xi_psi'), zlib=True)
    v.long_name = "mask on PSI-points" 
    v.option_0 = "land" 
    v.option_1 = "water" 

    # These variables, xl and el are not used by ROMS
    # but are required?

    v = nc.createVariable('xl', 'd', ())
    v.long_name = "domain length in the XI-direction" 
    v.units = "meter" 

    v = nc.createVariable('el', 'd', ())
    v.long_name = "domain length in the ETA-direction" 
    v.units = "meter" 

    nc.close()


def create_grid(gmap, grid_name, file_name='', 
                global_attributes={},
                format='NETCDF3_CLASSIC'):
    """
    Create a new ROMS grid file for a polar stereographic grid

    Arguments:
      gmap      : a gridmap.PolarStereographic instance
      grid_name : name of the grid
      file_name : name of the grid file,
                  default = '' giving grid_name + '_grid.nc'
      format    : 'NETCDF3_CLASSIC' or 'NETCDF4_CLASSIC'
                  default = 'NETCDF3_CLASSIC'

    Fills in geometric variables (lons, lats, metric, Coriolis).
    Makes space for topographic variables (h, hraw, masks).
    Also makes coordinate variables and includes grid mapping info
    following the CF-standard.
    
    """

    if not file_name:  # Use default
        file_name = grid_name + '_grid.nc'

    make_empty_gridfile(grid_name, file_name, gmap.Lm, gmap.Mm, 
                        global_attributes=global_attributes,
                        format=format)

    nc = Dataset(file_name, 'a')

    gridmap_varname = 'grid_mapping' # Name of grid mapping variable
    Lm, Mm = gmap.Lm, gmap.Mm

    # --- Grid map

    #v = nc.createVariable(gridmap_varname, 'i', ())
    v = nc.variables[gridmap_varname]
    #v.long_name = "grid mapping"
    d = gmap.CFmapping_dict()
    for att in d:
        setattr(v, att, d[att])
    v.proj4string = gmap.proj4string


    # ------------------------------------------------------
    # Compute variables defined by only by the grid mapping
    # ------------------------------------------------------

    #print "Saving geometric variables"

    # -----------------------
    # Coordinate variables
    # -----------------------
    nc.variables['xi_rho'][:]  = gmap.dx*np.arange(Lm+2)
    nc.variables['eta_rho'][:] = gmap.dx*np.arange(Mm+2)
    nc.variables['xi_u'][:]    = gmap.dx*(np.arange(Lm+1)+0.5)
    nc.variables['eta_u'][:]   = gmap.dx*np.arange(Mm+2)
    nc.variables['xi_v'][:]    = gmap.dx*np.arange(Lm+2)
    nc.variables['eta_v'][:]   = gmap.dx*(np.arange(Mm+1)+0.5)

    # ----------
    # Vertices 
    # ----------

    # Vertices at every half point in the grid
    # -0.5, 0, 0.5, ...., Lm+1.5
    Lvert = 2*Lm + 5
    Mvert = 2*Mm + 5
    X0 = 0.5*np.arange(Lvert)-0.5
    Y0 = 0.5*np.arange(Mvert)-0.5
    # Make 2D arrays with grid coordonates
    Xvert, Yvert = np.meshgrid(X0, Y0)
    Xrho = Xvert[1::2, 1::2]
    Yrho = Yvert[1::2, 1::2]

    lon_vert, lat_vert = gmap.grid2ll(Xvert, Yvert)

    # Set the different points
    nc.variables['lon_rho'][:,:] = lon_vert[1::2, 1::2]
    lat_rho = lat_vert[1::2, 1::2]
    nc.variables['lat_rho'][:,:] = lat_rho
    nc.variables['lon_u'][:,:]   = lon_vert[1::2, 2:-1:2]
    nc.variables['lat_u'][:,:]   = lat_vert[1::2, 2:-1:2]
    nc.variables['lon_v'][:,:]   = lon_vert[2:-1:2, 1::2]
    nc.variables['lat_v'][:,:]   = lat_vert[2:-1:2, 1::2]

    # ----------------------
    # Metric coefficients
    # ----------------------
    
    #pm = 1.0 / (gmap.map_scale(Xrho, Yrho) * gmap.dx)
    pm = gmap.map_scale(Xrho, Yrho) / gmap.dx
    pn = pm
    nc.variables['pm'][:,:] = pm
    nc.variables['pn'][:,:] = pn

    # Alternative: 
    # Could define pm and pn by differencing on the ellipsoid
    # However, for WGS84 the distance formula is complicated

    # --- Derivatives of metric coefficients

    # Use differencing, as the calculus is complicated for WGS84
    # the pm and pn fields are changing slowly and no
    # problems near the North Pole

    dndx = np.zeros_like(pm)
    dmde = np.zeros_like(pm)

    # Central differences
    dndx[:, 1:-1] = 0.5/pn[:, 2:] - 0.5/pn[:, :-2]
    dmde[1:-1, :] = 0.5/pm[2:, :] - 0.5/pm[:-2, :]

    # linear extrapolation to boundary
    dndx[:,0]  = 2*dndx[:,1]  - dndx[:,2]
    dndx[:,-1] = 2*dndx[:,-2] - dndx[:,-3]
    dmde[0,:]  = 2*dmde[1,:]  - dmde[2,:]
    dmde[-1,:] = 2*dmde[-2,:] - dmde[-3,:]

    # Alternative for spherical earth
    #phi0 = gmap.lat_ts*np.pi/180.0
    #R = gmap.ellipsoid.a
    #dndx = - (Xrho - gmap.xp)*gmap.dx / (pn**2 * R**2 * (1+np.sin(phi0)))
    #dmde = - (Yrho - gmap.yp)*gmap.dx / (pm**2 * R**2 * (1+np.sin(phi0)))

    # save the coefficients
    nc.variables['dndx'][:,:] = dndx
    nc.variables['dmde'][:,:] = dmde

    # ---------
    # Coriolis 
    # ---------

    Aomega = 2 * np.pi * (1+1/365.24) / 86400 # earth rotation
    nc.variables['f'][:,:] = 2 * Aomega * np.sin(lat_rho*np.pi/180.0)

    # ----------------
    # Rotation angle
    # ----------------

    nc.variables['angle'][:,:] = gmap.angle(Xrho, Yrho)
    # Could also be computed by differencing,
    # this would be very inaccurate near the North Pole

    # ------------------
    # Misc. variables
    # ------------------

    nc.variables['spherical'].assignValue('T')
    nc.variables['xl'].assignValue((Lm+1)*gmap.dx)
    nc.variables['el'].assignValue((Mm+1)*gmap.dx)

    # ---------------------
    # Close the grid file
    # ---------------------

    nc.close()

# -------------------------------------------------------

def subgridfile(file0, file1, i0, j0, Lm, Mm):
    ### Funker ikke helt

    f0 = Dataset(file0)
    gmap0 = gridmap.fromfile(f0)

    gmap1 = gridmap.subgrid(gmap0, i0, j0, Lm, Mm)

    grid_name = f0.gridname + "_sub"

    gridmap_varname = "grid_mapping"  # Les denne fra filen

    # Make an empty grid file of the correct shape
    make_empty_gridfile(grid_name, file1, Lm, Mm, format=f0.file_format)

    # Open this grid file
    f1 = Dataset(file1, 'a')
    # Add grid mapping
    v = f1.variables[gridmap_varname]
    #v.long_name = "grid mapping"
    d = gmap1.CFmapping_dict()
    for att in d:
        setattr(v, att, d[att])
    v.proj4string = gmap1.proj4string

    # -----------------------
    # Coordinate variables
    # -----------------------
    f1.variables['xi_rho'][:]  = gmap1.dx*np.arange(Lm+2)
    f1.variables['eta_rho'][:] = gmap1.dx*np.arange(Mm+2)
    f1.variables['xi_u'][:]    = gmap1.dx*(np.arange(Lm+1)+0.5)
    f1.variables['eta_u'][:]   = gmap1.dx*np.arange(Mm+2)
    f1.variables['xi_v'][:]    = gmap1.dx*np.arange(Lm+2)
    f1.variables['eta_v'][:]   = gmap1.dx*(np.arange(Mm+1)+0.5)

    # -------------------------
    # rho-point variables
    # ------------------------

    vars = ['lon_rho', 'lat_rho', 'mask_rho', 
            'pm', 'pn', 'dmde', 'dndx', 'angle',
            'f', 'h']

    for var in vars:
        v0 = f0.variables[var]
        v1 = f1.variables[var]
        v1[:,:] = v0[j0:j0+Mm+2, i0:i0+Lm+2]

    # hraw is special
    v0 = f0.variables['hraw']
    v1 = f1.variables['hraw']
    for t in range(len(f0.dimensions['bath'])):
        v1[t,:,:] = v0[j0:j0+Mm+2, i0:i0+Lm+2]

    # u-point variables
    vars = ['lon_u', 'lat_u', 'mask_u']
    for var in vars:
        v0 = f0.variables[var]
        v1 = f1.variables[var]
        # Eller er det i0+1:i0+Lm+2
        v1[:,:] = v0[j0:j0+Mm+2, i0:i0+Lm+1]
        
    # v-point variables
    vars = ['lon_v', 'lat_v', 'mask_v']
    for var in vars:
        v0 = f0.variables[var]
        v1 = f1.variables[var]
        v1[:,:] = v0[j0:j0+Mm+1, i0:i0+Lm+2]

    # psi-point variables
    vars = ['mask_psi']
    for var in vars:
        v0 = f0.variables[var]
        v1 = f1.variables[var]
        # sjekk om indeksering er forskjøvet
        v1[:,:] = v0[j0:j0+Mm+1, i0:i0+Lm+1]

    # Some special variables
    f1.variables['spherical'].assignValue('T')
    f1.variables['xl'].assignValue((Lm+1)*gmap1.dx)
    f1.variables['el'].assignValue((Mm+1)*gmap1.dx)

    f1.close()



    
