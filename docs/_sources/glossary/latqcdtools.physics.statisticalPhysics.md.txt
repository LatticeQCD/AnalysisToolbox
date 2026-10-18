latqcdtools.physics.statisticalPhysics
=============

```python
_compareWithZero(err, tol) -> bool:
"""
Check that err is consistent with zero. If err has an error bar, require that the probability of
getting a deviation at least this large by chance is above 0.01, as in _statisticalConsistencyCheck.
Otherwise require |err| < tol.
"""
```
```python
_getParameter(parameterInfo):
"""
Give the weighted average of literature values of some parameter. These are
in units with k_B=J=1 and no magnetic field. If there is only one known value,
give that back.

We follow the PDG prescription for averaging measurements x_i +/- s_i, see Section 5.2.2
of the introduction to the Review of Particle Physics, DOI: 10.1103/PhysRevD.110.030001.
With weights w_i = 1/s_i^2, the average and its error are
    xbar = sum_i w_i x_i / sum_i w_i,    s = (sum_i w_i)^(-1/2).
If the x_i disagree, s underestimates the uncertainty. Therefore compute
    chi^2 = sum_i w_i (x_i - xbar)^2,    S = sqrt( chi^2/(N-1) ),
and if S > 1, scale s up to S*s. The PDG computes S using only the x_i with s_i < 3 sqrt(N) s;
we use all of them. The weights are plain numbers, so xbar is a linear combination of the
gvars x_i, and gvar still propagates errors into exponents derived from xbar.
"""
```
```python
_printExponent(prefix, exponent)
```
```python
_statisticalConsistencyCheck(parameterInfo, obs)
```
```python
reweight(X, pRW, p0, S):
"""
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
"""
```
```python
class O2_3d():
```
```python
class O3_3d():
```
```python
class O4_3d():
```
```python
class Oinf_3d():
```
```python
class S3_2d():
"""
Exact solution for 2d S_3 class from Baxter "Exactly Solved Models in Statistical Mechanics"
"""
```
```python
class S4_2d():
"""
Exact solution for 2d S_4 class from Baxter, "Exactly Solved Models in Statistical Mechanics"
"""
```
```python
class UniversalityClass():
"""
Skeleton universality class from which all others inherit.
"""
```
```python
class Z2_2d():
"""
Onsager solution for 2d Z_2 class. 
"""
```
```python
class Z2_3d():
```
