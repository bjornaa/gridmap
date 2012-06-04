# -*- coding: utf-8 -*

"""Defines grid mappings, between longitude, latitude and grid coordinates"""

from numpy import pi, sqrt, sin, cos, tan, arctan, arctan2
try:
    from mpl_toolkits.basemap import Basemap
    has_basemap = True
except:
    has_basemap = False

all = ['Ellipsoid', 'Sphere', 'sphere', 'WGS84', 
       'PolarStereographic', 'gridmap_fromfile']

deg = 180.0 / pi
rad = pi / 180.0

# -------------------------
# Ellipsoid classes
# -------------------------

class Ellipsoid(object):
    def __init__(self, a, invf=None, name='general'):
        self.a = a                        # major semi-axis [m]
        self.invf = invf                  # inverse flattening
        if invf == None:    # sphere case
            self.f = 0.0                  # flatening  
        else:
            self.f = 1.0 / invf           
        self.e = sqrt((2-self.f)*self.f)  # eccentrity
        self.b = a*(1-self.f)             # minor semi-axis [m]
        self.name = name

class Sphere(Ellipsoid):
    def __init__(self, radius=6371000.0):  # met.no default
        Ellipsoid.__init__(self, a=radius, name='Sphere')
                
# Common instances
sphere = Sphere()   # default radius
WGS84  = Ellipsoid(a=6378137.0, invf=298.257223563, name='WGS84')

# -------------------------------------
# Polar Stereographic grid map class
# -------------------------------------

class PolarStereographic(object):
    """Polar stereographic grid mapping"""

    def __init__(self, xp, yp, dx, ylon,
                 shape=None, lat_ts=None, ellipsoid=sphere):
        self.xp    = xp     # x-coordinate of north pole
        self.yp    = yp     # y-coordinate of north pole
        self.dx    = dx     # grid resolution [m]
        self.ylon  = ylon   # longitude parallell to y-axis [deg]

        self.shape = shape  # Number of internal grid cells (Lm, Mm)
        
        if lat_ts:
            self.lat_ts = lat_ts  # latitude of true scale [deg]
        else:
            self.lat_ts = 60.0   # met.no default 
        self.ellipsoid = ellipsoid

        phi_c = self.lat_ts*rad    
        e = self.ellipsoid.e
        t_c = tan(0.25*pi-0.5*phi_c)             \
               / ((1-e*sin(phi_c))/(1+e*sin(phi_c)))**(0.5*e)
        m_c = cos(phi_c) / sqrt(1-(e*sin(phi_c))**2)
        #self.k_0 = 0.5 * m_c * sqrt((1+e)**(1+e)*(1-e)**(1-e)) / t_c
        self.t_c, self.m_c = t_c, m_c

        # Make an option string for proj4
        self.projstring = self._projstring()

        
    def ll2grid(self, lon, lat):
        """Compute grid coordinates x, y from longitude and latitude"""
        lambda0 = self.ylon*rad
        lambda_ = lon*rad
        phi     = lat*rad
        sinphi  = sin(phi)
        e = self.ellipsoid.e
        a = self.ellipsoid.a

        t = tan(0.25*pi-0.5*phi) / ((1-e*sinphi)/(1+e*sinphi))**(0.5*e)
        rho = a * t * self.m_c / (self.dx * self.t_c)
        x = self.xp + rho*sin(lambda_ - lambda0)
        y = self.yp - rho*cos(lambda_ - lambda0)
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

        rho = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp))
        t = rho * self.t_c * self.dx / (a * self.m_c)
        chi = 0.5*pi - 2*arctan(t)   # Conformal latitude
        if e > 0:
            phi = chi                                                         \
                 + (e**2/2 + 5*e**4/24 + e**6/12 + 13*e**8/360) * sin(2*chi)  \
                 + (7*e**4/48 + 29*e**6/240 + 811*e**8/11520) * sin(4*chi)    \
*                 + (7*e**6/120 +  81*e**8/1120) * sin(6*chi)                  \
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
        rho = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp)) * self.dx
        m = cos(phi) / sqrt(1-(e*sin(phi))**2)
        return rho / (a * m)

    def _projstring(self):
        proj = "+proj=stere"
        ellipsoid = self.ellipsoid
        if ellipsoid.name == 'WGS84':
            ellps = "+ellps=WGS84" 
        else:
            ellps = "+R=" + str(self.ellipsoid.a)
        lat_0 = "+lat_0=90"
        lon_0 = "+lon_0=" + str(self.ylon)
        x_0 = "+x_0=" + str(self.xp*self.dx)
        y_0 = "+y_0=" + str(self.yp*self.dx)
        lat_ts = "+lat_ts=" + str(self.lat_ts)
        projlist = [proj, ellps, lat_0, lat_ts, x_0, y_0, lon_0]
        return " ".join(projlist)

    
    if has_basemap:
        def basemap(self, Lm, Mm, resolution='i', area_thresh=None):
            # Basere det hele på 0 til Lm+1
            # => projeksjon = samme, plott litt begrenset
            #                bare halve randceller

            # Check ellipsoid
            ellipsoid = self.ellipsoid
            if ellipsoid.name == 'WGS84':
                rsphere = (ellipsoid.a, ellipsoid.b)
            else:
                rsphere = ellipsoid.a

            if area_thresh == None:
                # Default = half grid cell area
                area_thresh = 0.5*(self.dx**2) / 1.e6  

            # Compute lon/lat of grid corners (needed by basemap)
            lon0, lat0 = self.grid2ll(-0.5, -0.5)
            lon1, lat1 = self.grid2ll(Lm+1.5, Mm+1.5)

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
    grid_mapping_variable = nc.variables[var].mapping
    v = nc.variables[grid_mapping_variable]
    # 
    if v.ellipsoid == 'WGS84':
        ellipsoid = WGS84
    else:
        ellipsoid = Sphere(v.earth_radius)
    xp = v.false_easting
    yp = v.false_northing
    ylon = v.straight_vertical_longitude_from_pole
    dx = v.dx
    lat_ts = v.standard_parallel
    Lp, Mp = nc.variables[var].shape
    shape = (Lp-2, Mp-2)

    return PolarStereographic(xp, yp, dx, ylon, shape lat_ts, ellipsoid)




    
    
   
             
    

