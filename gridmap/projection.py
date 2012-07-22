# -*- coding: utf-8 -*

"""Module for grid map projections for ROMS grid
"""

# --------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# 2012-07-20
# --------------------------------

from numpy import pi, sqrt, sin, cos, tan, arctan, arctan2
try:
    from mpl_toolkits.basemap import Basemap
    has_basemap = True
except:
    has_basemap = False

__all__ = ['Ellipsoid', 'Sphere', 'sphere', 'WGS84', 
       'PolarStereographic', 'fromfile', 'rotate', 'subgrid']

deg = 180.0 / pi
rad = pi / 180.0

# -------------------------
# Ellipsoid classes
# -------------------------

class Ellipsoid(object):

    """Earth ellipsoid"""

    def __init__(self, a, invf=None):
        """Ellipsoid parameters

           a : major semi-axis [m]
           invf : optional, inverse flattening
                  defailt None for sphere

        """
        self.a = a                        # major semi-axis [m]
        self.invf = invf                  # inverse flattening
        if invf == None:    # sphere case
            self.f = 0.0                  # flatening  
        else:
            self.f = 1.0 / invf           
        self.e = sqrt((2-self.f)*self.f)  # eccentrity
        self.b = a*(1-self.f)             # minor semi-axis [m]

class Sphere(Ellipsoid):
    def __init__(self, radius=6371000.0):  # met.no default
        Ellipsoid.__init__(self, a=radius)
                
# Common instances
sphere = Sphere()   # default radius
WGS84  = Ellipsoid(a=6378137.0, invf=298.257223563)

# -------------------------------------
# Polar Stereographic grid map class
# -------------------------------------

class PolarStereographic(object):
    """Polar stereographic grid mapping"""

    def __init__(self, xp, yp, dx, ylon,
                 Lm=0, Mm=0, ellipsoid="sphere", lat_ts=60.0):
        self.xp    = xp     # x-coordinate of north pole
        self.yp    = yp     # y-coordinate of north pole
        self.dx    = dx     # grid resolution [m]
        self.ylon  = ylon   # longitude parallell to y-axis [deg]
        self.Lm = Lm        # Number of internal grid points, x-direction
        self.Mm = Mm        # Number of internal grid points, y-direction
        self.lat_ts = lat_ts  # latitude of true scale [deg]
        # Allow string arguments for ellipsoid
        if ellipsoid == "WGS84":  ellipsoid = WGS84
        if ellipsoid == "sphere": ellipsoid = sphere
        self.ellipsoid = ellipsoid

        phi_c = self.lat_ts*rad    
        e = self.ellipsoid.e
        self.t_c = tan(0.25*pi-0.5*phi_c)             \
               / ((1-e*sin(phi_c))/(1+e*sin(phi_c)))**(0.5*e)
        self.m_c = cos(phi_c) / sqrt(1-(e*sin(phi_c))**2)

        # Make an option string for proj4
        self.proj4string = self._proj4string()

        
    def ll2grid(self, lon, lat):
        """Compute grid coordinates x, y from longitude and latitude"""
        lambda0 = self.ylon*rad
        lambda_ = lon*rad
        phi     = lat*rad
        sinphi  = sin(phi)
        e = self.ellipsoid.e
        a = self.ellipsoid.a

        t = tan(0.25*pi-0.5*phi) / ((1-e*sinphi)/(1+e*sinphi))**(0.5*e)
        r = a * t * self.m_c / (self.dx * self.t_c)
        x = self.xp + r*sin(lambda_ - lambda0)
        y = self.yp - r*cos(lambda_ - lambda0)
        return x, y

    def grid2ll(self, x, y):
        """Compute longitude and latitude from grid coordinates x and y"""
        lambda0 = self.ylon*rad
        xp, yp = self.xp, self.yp
        a = self.ellipsoid.a
        e = self.ellipsoid.e
        
        lambda_= lambda0 + arctan2(x - xp, yp - y)
        # normalize to (-180,180]
        lambda_ = pi - (pi-lambda_) % (2*pi)

        r = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp))
        t = r * self.t_c * self.dx / (a * self.m_c)
        chi = 0.5*pi - 2*arctan(t)   # Conformal latitude
        if e > 0:
            phi = chi                                                         \
                 + (e**2/2 + 5*e**4/24 + e**6/12 + 13*e**8/360) * sin(2*chi)  \
                 + (7*e**4/48 + 29*e**6/240 + 811*e**8/11520) * sin(4*chi)    \
                 + (7*e**6/120 +  81*e**8/1120) * sin(6*chi)                  \
                 + (4279*e**8/161280) * sin(8*chi)
        else:
            phi = chi
        return (lambda_*deg, phi*deg)

    def angle(self, x, y):
        """Angle (in radians) from local y-direction anticlockwise to north"""
        
        return -arctan2(x-self.xp, self.yp-y)

    def map_scale(self, x, y):
        """Projection map scaling factor"""
        
        a = self.ellipsoid.a
        e = self.ellipsoid.e
        xp, yp = self.xp, self.yp

        lon, lat = self.grid2ll(x, y)
        phi = lat*rad
        r = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp)) * self.dx
        m = cos(phi) / sqrt(1-(e*sin(phi))**2)
        return r / (a * m)

    def CFmapping_dict(self):
        """Make a NetCDF grid map attribute dictionary

            Follows the CF standard, with extra attributes
            `ellipsoid` and `dx`
        """
        
        d = dict(grid_mapping_name = 'polar_stereographic',
                 latitude_of_projection_origin = 90.0,
                 straight_vertical_longitude_from_pole = self.ylon,
                 standard_parallel = self.lat_ts,
                 false_easting = self.xp*self.dx,
                 false_northing = self.yp*self.dx,
                 dx = self.dx)
                     
        if self.ellipsoid.invf: # Not sphere => WGS84
            d['ellipsoid'] = 'WGS84'
            d['semi_major_axis'] = self.ellipsoid.a
            d['inverse_flattening'] = self.ellipsoid.invf
        else: # ellipsoid == 'sphere'
            d['ellipsoid'] = 'sphere'
            d['earth_radius'] = self.ellipsoid.a

        return d
        

    def _proj4string(self):
        """Make an option string for proj4"""

        # proj4 representation of the ellipsoid
        if self.ellipsoid.invf == None:  # sphere
            ellps = "+R=" + str(self.ellipsoid.a)
        else:  # only other case is WGS84
            ellps = "+ellps=WGS84" 

        # The rest of the options are found in the CF-standard
        templ = ('+proj=stere',
                 ellps,
                 '+lat_0=%(latitude_of_projection_origin)s',
                 '+lat_ts=%(standard_parallel)s',
                 '+x_0=%(false_easting)s',
                 '+y_0=%(false_northing)s',
                 '+lon_0=%(straight_vertical_longitude_from_pole)s')
        template = ' '.join(templ)
        return template % self.CFmapping_dict()

    if has_basemap:
        def basemap(self, resolution='i', area_thresh=None):

            # Check ellipsoid
            ellipsoid = self.ellipsoid
            if ellipsoid.invf == None:  # sphere case
                rsphere = ellipsoid.a
            else:                       # ellipsoid case
                rsphere = (ellipsoid.a, ellipsoid.b)

            if area_thresh == None:
                # Default = half grid cell area
                area_thresh = 0.5*(self.dx**2) / 1.e6  

            # Compute lon/lat of grid corners (needed by basemap)
            lon0, lat0 = self.grid2ll(0.0, 0.0)
            Lm, Mm = self.Lm, self.Mm
            lon1, lat1 = self.grid2ll(Lm+1.0, Mm+1.0)
                
            return Basemap(llcrnrlon=lon0, llcrnrlat=lat0,
                           urcrnrlon=lon1, urcrnrlat=lat1,
                           projection='stere',           
                           resolution=resolution,
                           area_thresh=area_thresh, 
                           rsphere=rsphere,
                           lat_0 = 90.0, lon_0 = self.ylon,
                           lat_ts = self.lat_ts)
        
# -------------

def fromfile(nc, var='h'):
    """Create a PolarStereographic instance from a netCDF file

    arguments:
    nc : open netcdf-file (netCDF4.Dataset instance)
    var, optional : variable in the netcdf file with mapping attribvute
               default is topography 'h'
    
    returns: PolarStereographic instance

    The mapping attributes must follow the CF-standard
    + additional attributes, ellipsoid and dx.
    Grid files created by define_grid are OK.

    Typical usage: 
    >>> nc = NetCDF4.Dataset(roms_grid_file)
    >>> gmap = gridmap_fromfile(nc)

    """

    # Legge til fornuftig feilkontroll
    #print "nc = ", nc
    #print "var = ", var
    #print "nc.variables[var] = ", nc.variables[var]
    grid_mapping_variable = nc.variables[var].grid_mapping
    v = nc.variables[grid_mapping_variable]
    # 
    if v.ellipsoid == 'WGS84':
        ellipsoid = WGS84
    else:
        ellipsoid = Sphere(v.earth_radius)
    dx = v.dx
    xp = v.false_easting / dx
    yp = v.false_northing / dx
    ylon = v.straight_vertical_longitude_from_pole
    lat_ts = v.standard_parallel
    Mp, Lp = nc.variables[var].shape
    Lm, Mm = Lp-2, Mp-2

    #return PolarStereographic(xp, yp, dx, ylon, shape, lat_ts, ellipsoid)
    return PolarStereographic(xp, yp, dx, ylon,
                              Lm, Mm, ellipsoid, lat_ts)

# -----------------------------------

def rotate(gmap):
    """
    Rotate a grid mapping 90 degrees counterclockwise

    New grid coordinates (x1, y1) are given by
      x1 = Mm + 1 - y; y1 = x
    """

    xp = gmap.Mm + 1 - gmap.yp
    yp = gmap.xp
    ylon = gmap.ylon - 90
    if ylon <= -180: ylon = ylon + 360  # Normalize to (-180,180]
    Lm = gmap.Mm
    Mm = gmap.Lm

    return PolarStereographic(xp, yp, gmap.dx, ylon, Lm, Mm,
                              gmap.ellipsoid, gmap.lat_ts)

# -----------

def subgrid(gmap, i0, j0, Lm, Mm):
    """
    Return a subgrid with origin at (i0,j0)
    
    New grid coordinates (x1, y1) are given by
      x1 = x - i0,  y1 = y - j0

    Subgrid requirement, must have
      i0 + Lm <= gmap.Lm; j0 + Mm <= gmap.Mm
    """

    if i0 + Lm > gmap.Lm:
        raise ValueError, "i0 + Lm too large"
    if j0 + Mm > gmap.Mm:
        raise ValueError, "j0 + Mm too large"
    
    xp = gmap.xp - i0
    yp = gmap.yp - j0

    return PolarStereographic(xp, yp, gmap.dx, gmap.ylon, Lm, Mm,
                              gmap.ellipsoid, gmap.lat_ts)
    

