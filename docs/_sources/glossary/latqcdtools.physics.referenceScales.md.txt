latqcdtools.physics.referenceScales
=============

```python
_betaRangeWarn(beta, beta_range):
"""
Many of these ansätze a(beta) have coefficients that were determined by performing a fit within a certain
beta range. This warning flashes whenever you are using information outside of that range, where the ansatz
is less likely be to reliable.

Args:
    beta (float) or numpy array
    beta_range (array-like): min and max beta of range, in that order 
"""
```
```python
a_div_r1(beta, year):
"""
Get a/r_1(beta).

Args:
    beta (float)
    year (int/str): year that parameterization was determined 

Returns:
    float: a/r_1
"""
```
```python
a_times_fk(beta, year):
"""
Get a*f_k(beta).

Args:
    beta (float)
    year (int/str): year that parameterization was determined 

Returns:
    float: a*f_k 
"""
```
```python
a_times_ms_2014(beta)
```
```python
allton_type_ansatz(beta, c0, c2, d2)
```
```python
fit_2014Eos_eqB2(beta, c0, c2, d2)
```
```python
fit_tayloraLambda(beta, a, b, c)
```
```python
ignoreBetaRange():
"""
Turn off the beta range warnings. 
"""
```
```python
r0_div_a(beta, year):
"""
Get r0/a(beta).

Args:
    beta (float)
    year (int/str): year that parameterization was determined 

Returns:
    float: r0/a 
"""
```
```python
r1_times_ms_2014(beta)
```
```python
sqrtt0_div_a(beta):
"""
Get sqrt(t0/a)(beta)

Args:
    beta (float)

Returns:
    float: sqrt(t0/a) 
"""
```
```python
wuppertal_type_ansatz(beta, c1, c2, c3, c4)
```
