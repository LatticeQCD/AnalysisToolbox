latqcdtools.physics.statisticalPhysics
=============

`_printExponent(prefix, exponent)`


`reweight(X, pRW, p0, S)`
 
    Reweight an observable X computed at a simulation point p0 to a nearby
    simulation point pRW. We assume the action depends linearly on the simulation
    parameter, i.e. S' ~ p S

    By the way, if you are going to pair this with a jackknife to estimate e.g.
    where a response function is maximized, make sure you pick enough reweighting
    points to properly resolve where the maximum is.

    Args:
        X (np.array): Measurements to reweight. 
        pRW (float): Reweight to this target. 
        p0 (float): Simulation point.
        S (np.array): Measurements of the action (extensive) divided by parameter p. 
    
`O2_3d()`
 
    3d O(2) critical exponents from JHEP08 (2016) 036 
    
`O3_3d()`
 
    3d O(3) critical exponents from JHEP08 (2016) 036 
    
`O4_3d()`
 
    3d O(4) critical exponents from Nucl. Phys. B 675, 533-554 (2003). 
    
`S3_2d()`
 
    Exact solution for 2d S_3 class from Baxter "Exactly Solved Models in Statistical Mechanics"
    
`S4_2d()`
 
    Exact solution for 2d S_4 class from Baxter, "Exactly Solved Models in Statistical Mechanics"
    
`UniversalityClass()`
 
    Skeleton universality class from which all others inherit.
    
`Z2_2d()`
 
    Exact solution for 2d Z_2 class. 
    
`Z2_3d()`
 
    3d Z_2 critical exponents from J. Stat. Phys. 157. 869-914 (2014). 
    
