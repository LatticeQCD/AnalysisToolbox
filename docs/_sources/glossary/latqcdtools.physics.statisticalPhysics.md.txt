latqcdtools.physics.statisticalPhysics
=============

`_printExponent(prefix, exponent)`


`reweight(X, pRW, p0, S)`
 
    Reweight an observable X computed at a simulation point p0 to a nearby
    simulation point pRW. We assume the action depends linearly on the simulation
    parameter, i.e. S' ~ p S

    Args:
        X (np.array): Measurements to reweight. 
        pRW (float): Reweight to this target. 
        p0 (float): Simulation point.
        S (np.array): Measurements of the action (extensive) divided by parameter p. 
    
`O2_3d()`
 
    3d O(2) critical exponents from Phys. Lett. B 492, 219 (2000). 
    
`O3_3d()`
 
    3d O(3) critical exponents from https://en.wikipedia.org/wiki/Universality_class. 
    
`O4_3d()`
 
    3d O(4) critical exponents from Nucl. Phys. B 675, 533-554 (2003). 
    
`UniversalityClass()`
 
    Skeleton universality class from which all others inherit.
    
`Z2_2d()`
 
    Exact solution for 2d Z_2 class. 
    
`Z2_3d()`
 
    3d Z_2 critical exponents from J. Stat. Phys. 157. 869-914 (2014). 
    
`Z2_4d()`
 
    Exact solution for 2d Z_2 class. 
    
