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

.. math:: dr = - \frac{r_0}{1 + \sin \phi} d\phi = 
               - \frac{r}{\cos \phi}  d\phi

Introduce the map factor :math:`m` by 

.. math:: dr = -m R d\phi

or explicitly,

.. math:: m = \frac{r_0}{R(1+\sin \phi)} 
            = \frac{r}{R \cos \phi}

The scale factor :math:`r_0` can be determined by requiring true scale
at latitude :math:`\phi_0`, i.e. :math:`m(\phi_0) = 1` or

.. math:: r_0 = R (1 + \sin \phi_0)

This throws out the earth radius, giving the expression

.. math:: m = \frac{1 + \sin \phi_0}{1 + \sin \phi}


The standard spherical metric is

.. math:: g &= R^2 \cos^2 \phi d\lambda \otimes d\lambda
             + R^2 d\phi \otimes d\phi                     \\
            &= R^2 \cos^2 \phi d\lambda \otimes d\lambda
             + \frac{R^2}{r_0^2} (1 + \sin \phi)^2 dr \otimes dr \\
            &= \frac{1}{m^2} 
              (r^2 d\lambda \otimes d\lambda + dr \otimes dr)    \\
            &= \frac{1}{m^2}(dx \otimes dx + dy \otimes dy)


The lack of cross term shows that this is an `orthogonal` coordinate
system and the same metric coefficient for both terms shows that it is
`conformal` (angle-preserving). ROMS needs the derivatives of the inverse
map factor.

.. math:: dm = \frac{1}{r_0 R} rdr = \frac{1}{r_0 R}(xdx + ydy)

.. math:: d(\frac{1}{m}) = - \frac{1}{m^2} dm
          = -\frac{1}{r_0 R m^2} (xdx + ydy)


.. math::  \frac{\partial}{\partial x} \left( \frac{1}{m} \right)
                   =  - \frac{x}{r_0 R m^2}, \quad
           \frac{\partial}{\partial y} \left( \frac{1}{m} \right)
                   =  - \frac{y}{r_0 R m^2}, \quad

In python code::

  import numpy as np 
  ...
  r = np.sqrt((i-xp)**2 + (j-yp)**2)
  r0 = R * (1+np.sin(phi0))
  ...
  pm[j,i] = r / (R * np.cos(phi[j,i]) * dx)
  pn[j,i] = pm[j,i]
  ...
  dndx[j,i] = (xp-i)*dx / (r0 * R * pn[j,i]**2)
  dmde[j,i] = (yp-j)*dx / (r0 * R * pm[j,i]**2)



The general ellipsoid case
--------------------------

The earth is modelled as an oblate rotational ellipsoid with
equatorial radius (major semiaxis) :math:`a` and  polar radius
(minor semiaxis) :math:`b`. The flattening :math:`f` and eccentricity
:math:`e` are defined by:

.. math:: f = 1 - \frac{b}{a}, \quad 
          e^2 = 1 - \frac{b^2}{a^2} = 2f - f^2

The two ellipsoids commonly used are the *sphere* with default radius
:math:`a = b = R = 6371000` m, or the *WGS84* ellipsoid defined by
:math:`a=6378137` m and :math:`f=1/298.257223563`.

The geometric projection from the South Pole to a plane intersecting
the ellipsoid at constant latititude :math:`\phi_0` will work, but is
not a *conformal* mapping. Instead the stereographic projection is
defined by

.. math:: r = r_0 \frac{\cos \phi}{1+\sin \phi}E(\phi)
            = r_0 \frac{\cos \chi}{1+\sin \chi}

where :math:`E` is a correction term

.. math:: E = \left( \frac{1-e\sin \phi}
                          {1+e\sin \phi} \right)^{\frac{e}{2}} 

and :math:`\chi = \chi(\phi)` is the *conformal* latitude.

.. math:: \chi = \frac{\pi}{2} 
      - 2 \arctan \left( \frac{E(\phi) \cos \phi}{1 + \sin \phi} \right )

The length element along a meridian is

.. math:: dS = \frac{a(1-e^2)}{(1-e^2\sin^2 \phi)^{3/2}} d\phi

and the map conversion factor is defined by

.. math:: dr = m dS, \quad \text{or} \quad m = ...

The requirement of true scale at :math:`\phi_0` gives

.. math:: r_0 = 



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

To obtain true scale at latitude :math:`\phi_c`, the scale factor
:math:`r_0` is given by

.. math:: r_0 =  a \frac{1 + \sin \phi_0}{\sqrt{1 - e^2 \sin^2 \phi_0}}
                   E(\phi_0)^{-1}
  
--- gammel tekst ---
   
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





