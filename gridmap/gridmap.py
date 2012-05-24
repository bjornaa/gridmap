# -*- coding: utf-8 -*

"""Defines grid mappings, between longitude, latitude and grid coordinates"""

from numpy import pi, sqrt, sin, cos, tan, arctan, arctan2

all = ['Ellipsoid', 'Sphere', 'sphere', 'WGS84', 
       'PolarStereographic', 'gridmap_fromfile']

deg = 180.0 / pi
rad = pi / 180.0

# -------------------------
# Ellipsoid classes
# -------------------------

class Ellipsoid(object):
    def __init__(self, a, invf=None):
        self.a = a                        # major semi-axis [m]
        self.invf = invf                  # inverse flattening
        if invf == None:    # sphere case
            self.f = 0.0                  # flatening  
        else:
            self.f = 1.0 / invf           
        self.e = sqrt((2-self.f)*self.f)  # eccentrity
        self.b = a*(1-self.f)             # minor semi-axis [m]
        #self.a = 6370997.0  # Med denne klaffer det
        # mer og mer feil jo lengre unna.

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

    def __init__(self, xp, yp, dx, ylon, lat_c=None, ellipsoid=sphere):
        self.xp    = xp     # x-coordinate of north pole
        self.yp    = yp     # y-coordinate of north pole
        self.dx    = dx     # grid resolution [m]
        self.ylon  = ylon   # longitude parallell to y-axis [deg]
        if lat_c:
            self.lat_c = lat_c  # latitude of true scale [deg]
        else:
            self.lat_c = 60.0   # met.no default 
        self.ellipsoid = ellipsoid

        phi_c = self.lat_c*rad    
        e = self.ellipsoid.e
        #print "e = ", e
        t_c = tan(0.25*pi-0.5*phi_c)             \
               / ((1-e*sin(phi_c))/(1+e*sin(phi_c)))**(0.5*e)
        m_c = cos(phi_c) / sqrt(1-(e*sin(phi_c))**2)
        self.k_0 = 0.5 * m_c * sqrt((1+e)**(1+e)*(1-e)**(1-e)) / t_c
        self.t_c, self.m_c = t_c, m_c
        #print "k_0 = ", self.k_0
        #print 0.5*cos(phi_c) / tan(0.25*pi - 0.5*phi_c)
        #print 0.5*(1+sin(phi_c))


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
        rho = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp)) * self.dx
        m = cos(phi) / sqrt(1-(e*sin(phi))**2)
        return rho / (a * m)

    def projstring1(self):
        """Return an option string for proj4"""
        projection = "+proj=stere"
        if isinstance(self.ellipsoid, Sphere):
            #ellipsis = "+ellps=sphere +a=" + str(self.ellipsoid.a)
            ellipsis = "+a=" + str(self.ellipsoid.a)
        else:  
            ellipsis = "+ellps=WGS84" 
        lat_0 = "+lat_0=90"
        lon_0 = "+lon_0=" + str(self.ylon)
        x_0 = "+x_0=" + str(self.xp)
        y_0 = "+y_0=" + str(self.yp)
        k_0 = "+k_0=" + str(self.k_0 / self.dx)
        
        projlist = [projection, ellipsis, lat_0, k_0,
            x_0, y_0, lon_0]

        return " ".join(projlist)

    def projstring2(self):
        """Return an alternative option string for proj4"""
        projection = "+proj=stere"
        if isinstance(self.ellipsoid, Sphere):
            #ellipsis = "+ellps=sphere +a=" + str(self.ellipsoid.a)
            ellipsis = "+a=" + str(self.ellipsoid.a)
        else:  
            ellipsis = "+ellps=WGS84" 
        lat_0 = "+lat_0=90"
        lon_0 = "+lon_0=" + str(self.ylon)
        x_0 = "+x_0=" + str(self.xp*self.dx)
        y_0 = "+y_0=" + str(self.yp*self.dx)
        lat_ts = "+lat_ts=" + str(self.lat_c)
        
        projlist = [projection, ellipsis, lat_0, lat_ts,
            x_0, y_0, lon_0]

        return " ".join(projlist)

# -------------

def gridmap_fromfile(nc, var='h'):
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
    lat_c = v.standard_parallel

    return PolarStereographic(xp, yp, dx, ylon, lat_c, ellipsoid)




    
    
   
             
    

