latqcdtools.physics.correlators
=============

```Python
_periodicLHS(corr, m, it):
'''
Helper function for finding effective mass curve, baking in the fact
that the lattice is periodic. See e.g. Gattringer and Lang eq. (6.57).
Meant to be used with a solver.
'''
```
```Python
effectiveMass(corr, algorithm='simple', guess=None) -> numpy.ndarray:
'''
Get m_eff(t). Where this plateaus can be used as a first estimate for
the ground state. See e.g. Gattringer and Lang eq. (6.56)

Args:
    corr (np.ndarray)
    algorithm (str): 'simple' or 'periodic'
    guess (float): initial guess for periodic mass meff[0]

Returns:
    np.ndarray: m_eff(t) 
'''
```
```Python
foldCorrelator(corr) -> numpy.ndarray:
'''
On a periodic lattice, information for t>Nt/2 is redundant. This method
combines corr below and above that threshold to improve statistics.

Args:
    corr (np.ndarray)

Returns:
    np.ndarray: folded correlator 
'''
```
