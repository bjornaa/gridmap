==================================
Polar stereographic map projection
==================================

The sphere case
---------------

Consider the earth as a sphere of radius :math:`R` and let
:math:`\phi` be the latitude and :math:`\lambda` the longitude.  The
Actic polar stereographic mapping with true latitude :math:`\phi_0` is
the projection from the South Pole to the plane intersecting the
sphere at this latitude, see figure X. The radius :math:`r` from the
North Pole to the projected point is given by:

.. math:: r(\phi) = r_0 \frac{\cos \phi}{1 + \sin \phi}
             = r_0 \tan \left( \frac{\pi}{4} - \frac{\phi}{2} \right)

Taking differentials gives

.. math:: dr = - \frac{R}{r_0}{1 + \sin \phi} d\phi

The standard spherical metric is

.. math:: g = R^2 \cos^2 \phi d\lambda \otimes d\lambda
             + R^2 d\phi \otimes d\phi

In the projected plane this becomes

.. math:: g = R^2 \cos^2 \phi d\lambda \otimes d\lambda
             + \frac{R^2}{r_0^2} (1 + \sin \phi)^2 dr \otimes dr

Using Cartesian coordinates :math:`x`, :math:`y` in the projected
plane, we have differential relations for polar coordinates

.. math:: dr = \frac{x}{r} dx + \frac{y}{r} dy, \\
          d\lambda = - \frac{y}{r} dx + \frac{x}{r} dy

This gives

.. math:: g = \frac{1}{m^2} dx \otimes dx + \frac{1}{n^2} dy \otimes dy

with metric coeffients

.. math:: m = n = \frac{r_0}{R (1 + \sin \phi)}

The lack of cross term shows that this is an `orthogonal` coordinate
system, and :math:`m=n` shows that it is `conformal`
(angle-preserving).  The scale factor :math:`r_0` can be determined by
requiring true scale at latitude :math:`\phi_0`, i.e. :math:`m(\phi_0) = 1` or

.. math:: r_0 = R (1 + \sin \phi_0)

This throws out the earth radius, giving the expression

.. math:: m = \frac{1 + \sin \phi_0}{1 + \sin \phi}

For computations in the orhogonal coordinates systems, the partial
derivatives :math:`m_y` and :math:`n_x` are needed. They are the
terms ``dmde`` and ``dndx`` in ROMS.

.. math:: d\frac{1}{m} = 
      - \frac{R}{r_0} \frac{\cos^2 \phi}{(1+\sin \phi)^2} dr
      = \frac{R}{r_0 m^2} (xdx + ydy)
            
(mellomledd er feil) This gives 

.. math:: \frac{\partial m^{-1}}{\partial x} = \frac{R x}{r_0 m^2} 
          = \frac{x(1 + \sin \phi)^2}{(1 + \sin \phi_0)^3}

and

.. math:: \frac{\partial n^{-1}}{\partial x} = 
          \frac{\partial m^{-1}}{\partial x} = \frac{R y}{r_0 m^2}


The general ellipsoid case
--------------------------

The earth is modelled as an oblate rotational ellipsoid with
equatorial radius (major semiaxis) :math:`a` and  polar radius
(minor semiaxis) :math:`b`. The flattening :math:`f` and eccentricity
:math:`e` are defined by:

.. math:: f = 1 - \frac{b}{a}, \quad e = 1 - \frac{b^2}{a^2} = 2f - f^2

Two ellipsoids are available, the *sphere* with default radius
:math:`a = b = R = 6371000` m, or the *WGS84* ellipsoid defined by
:math:`a=6378137` m and :math:`f=1/298.257223563`.


The projection is an *azimuthal* projection best described with polar
coodinates.

.. math:: \lambda, \phi \mapsto \lambda, \rho(\phi)

where :math:`\lambda` is longitude, :math:`\phi` is latitude and 
:math:`\rho` is the radius in the target plane. The polar
stereographic projection is the only *conformal* (angle-preserving)
projection of this kind. For an ellipsoid the projection is given by


.. math::  \rho(\phi) =  
              \rho_0 \tan \left( \frac{\pi}{4}-\frac{\phi}{2} \right) 
                  \left( \frac{1+e\sin \phi}
                            {1-e\sin \phi} \right)^{\frac{e}{2}} 

To obtain true scale at latitude :math:`\phi_c`, the scale factor :math:`\rho_0` is given by

.. math::  \rho_0 =  a \frac{\cos \phi_c}{\sqrt{1-e^2\sin^2 \phi_c}}
                    \tan \left( \frac{\pi}{4}+\frac{\phi_c}{2} \right)
                  \left( \frac{1-e\sin \phi_c}
                              {1+e\sin \phi_c} \right)^{\frac{e}{2}} 
   
For an ellipsoid this is impossible to invert analytically. The
inverse can be found by iteration (as done in *proj4*) or by the
following approximate series development used in
*gridmap.PolarStereographic.grid2ll*

.. math:: \phi \approx \chi + 
    \left( \frac{1}{2}e^2 + \frac{5}{24}e^4 + \frac{1}{12}e^6 +
            \frac{13}{360}e^8 \right) \sin 2\chi
    + \left( \frac{7}{48}e^4 + \frac{29}{240}e^6 + 
            \frac{811}{11520}e^8 \right) \sin 4\chi 
    + \left( \frac{7}{120}e^6 + \frac{81}{1120}e^8 \right) \sin 6\chi
    + \frac{4279}{161280}e^8  \sin 8\chi

where :math:`\chi` is the *conformal* latitude given by

.. math:: \rho = \rho_0 
            \tan \left( \frac{\pi}{4} - \frac{\chi}{2} \right), \quad
            \chi = \frac{\pi}{2} - 2 \arctan \frac{\rho}{\rho_0}
           

With :math:`\phi_c = 60^{\circ}\mathrm{N}` the error for station M with
longitude :math:`\phi = 66^{\circ}\mathrm{N}` the error is ...

For a sphere, the conformal latitude and the ordinary latitude are
equal and the projection simplifies to 

.. math:: \rho(\phi) =  \rho_0 
            \tan \left( \frac{\pi}{4} - \frac{\phi}{2} \right)
             = \rho_0 \frac{\cos \phi}{1 + \sin \phi}

with

.. math:: \rho_0 = a \cos \phi_c 
             \tan \left( \frac{\pi}{4} + \frac{\phi_c}{2} \right)
	     = a (1 + \sin \phi_c)  

and has a simple geometric interpretation as a projection from the
South pole onto the plane at the true latitude :math:`\phi_c`.

For more info on map projections, a standard reference is
Snyder(1987) 

Grid in the projected plane
---------------------------

A grid in the projected plane is given by xp, yp, dx, ylon where
xp and yp are the grid coordinates of the north pole, dx is the grid
spacing (at latitude of true scale) and ylon is the longitude
parallel to the y-axis. (see figure ...). Denote these quantitities 
by :math:`x_0, y_0, \Delta x, \lambda_0` in mathematical notation.

The grid coordinates are then computed by

.. math:: 
      x = x_0 + \frac{\rho(\phi)}{\Delta x} \sin( \lambda - \lambda_0)

      y = y_0 - \frac{\rho(\phi)}{\Delta x} \cos( \lambda - \lambda_0)


The map factor ..., angle variable




J.P. Snyder, 1987, Map Projections -- A Working Manual, 
US Geological Survey professional paper 1395
http://pubs.er.usgs.gov/publication/pp1395
direct link:
http://pubs.usgs.gov/pp/1395/report.pdf





