latqcdtools.physics.statisticalPhysics
=============

```Python
_compareWithZero(err, tol) -> bool
```
```Python
_getParameter(parameterInfo):
'''
Give the weighted average of literature values of some parameter. These are 
in units with k_B=J=1 and no magnetic field. If there is only one known value, 
give that back. 
'''
```
```Python
_printExponent(prefix, exponent)
```
```Python
_statisticalConsistencyCheck(parameterInfo, obs)
```
```Python
ignoreTension():
'''
Turn off statistical tension warnings. 
'''
```
```Python
reweight(X, pRW, p0, S):
'''
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
'''
```
```Python
class O2_3d():
```
```Python
class O3_3d():
```
```Python
class O4_3d():
```
```Python
class Oinf_3d():
```
```Python
class S3_2d():
'''
Exact solution for 2d S_3 class from Baxter "Exactly Solved Models in Statistical Mechanics"
'''
```
```Python
class S4_2d():
'''
Exact solution for 2d S_4 class from Baxter, "Exactly Solved Models in Statistical Mechanics"
'''
```
```Python
class UniversalityClass():
'''
Skeleton universality class from which all others inherit.
'''
```
```Python
class Z2_2d():
'''
Onsager solution for 2d Z_2 class. 
'''
```
```Python
class Z2_3d():
```
