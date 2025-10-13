latqcdtools.physics.continuumExtrap
=============

```Python
_powerSeries(x, coeffs):
'''
The default fit form for a continuum extrapolation is a power series in a^2.

Args:
    x (array-like): a^2 or 1/Nt^2 
    coeffs (array-like): power series coefficients 

Returns:
    array-like: power series in x 
'''
```
```Python
continuumExtrapolate(x, obs, obs_err, order=1, show_results=False, plot_results=False, prior=None, start_coeffs=None, priorsigma=None, error_strat='propagation', xtype='a', nproc=6, detailedInfo=False):
'''
A convenience wrapper for the Extrapolator. 
'''
```
```Python
class Extrapolator(x, obs, obs_err, order=1, xtype='a', error_strat='propagation', ansatz=None, nproc=6, tol=1e-12, max_fev=None):
```
