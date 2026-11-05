latqcdtools.physics.correlators
=============

```python
_periodicLHS(corr, m, it, Nt):
"""
Helper function for finding effective mass curve, baking in the fact
that the lattice is periodic. See e.g. Gattringer and Lang eq. (6.57).
Meant to be used with a solver.
"""
```
```python
effectiveMass(corr, algorithm='simple', guess=None, Nt=None) -> numpy.ndarray:
"""
Get m_eff(t). Where this plateaus can be used as a first estimate for
the ground state. See e.g. Gattringer and Lang eq. (6.56). The 'simple'
algorithm ignores the backward-propagating part of the correlator, so it
is biased low when exp(-m*(Nt-2t)) is not negligible; 'periodic' accounts for it.

Args:
    corr (np.ndarray)
    algorithm (str): 'simple' or 'periodic'
    guess (float): initial guess for periodic mass meff[0]
    Nt (int, optional): temporal extent of the lattice. Needed with 'periodic' when corr
      is folded, i.e. has Nt//2+1 entries. Defaults to len(corr), i.e. an unfolded corr.

Returns:
    np.ndarray: m_eff(t)
"""
```
```python
foldCorrelator(corr) -> numpy.ndarray:
"""
On a periodic lattice, information for t>Nt/2 is redundant. This method
combines corr below and above that threshold to improve statistics.

Args:
    corr (np.ndarray)

Returns:
    np.ndarray: folded correlator 
"""
```
