==================================
Polar stereographic map projection
==================================

Ellipsoid
---------

The earth is modelled as an oblate rotational ellipsoid with
equatorial radius (major semiaxis) :math:`a` and  polar radius
(minor semiaxis) :math:`b`. The flattening :math:`f` and eccentricity
:math:`e` are defined by:

.. math:: f = 1 - \frac{b}{a}, \quad e = 1 - \frac{b^2}{a^2}

Two ellipsoids are available, the *sphere* with default radius
:math:`a = b = 6371000` m, or the *WGS84* ellipsoid defined by
:math:`a=6378137` m and :math:`f=1/298.257223563`.

Projection
----------

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

To obtain true scale at latitude :math:`lon_0`, the scale factor :math:`\rho_0` is given by



For a sphere, the conformal latitude and the ordinary latitude are
equal and the projection simplifies to

.. math:: \rho(\phi) =  \tan \left( \frac{\pi}{4} - \frac{\phi}{2} \right)

and has a simple geometric interpretation as a projection from the
South pole onto the tangent plane at the North pole.


