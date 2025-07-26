latqcdtools.math.spline
=============

```Python
_even_knots(xdata, nknots):
'''
Return a list of nknots evenly spaced knots. 
'''
```
```Python
_random_knots(xdata, nknots, randomization_factor=1, SEED=None):
'''
Return a list of nknots randomly spaced knots. 
'''
```
```Python
getSpline(xdata, ydata, num_knots=None, edata=None, order=3, rand=False, fixedKnots=None, getAICc=False, natural=False):
'''
This is a wrapper that calls SciPy spline interpolation methods, depending on your needs. Generally
this uses scipy.interpolate.splrep, which uses B-splines. If natural=True and edata=None, it will
use scipy.interpolate.CubicSpline to solve. If natural=True and edata are provided, it will do a
smoothing spline that attempts to force no curvature at the endpoints, based on the penultimate points. 

Args:
    xdata (array-like)
    ydata (array-like)
    num_knots (int):
        The number of knots.
    edata (array-like, optional): 
        Error data. Defaults to None.
    order (int, optional):
        Order of the spline. Defaults to 3.
    rand (bool, optional): 
        Use randomly placed knots? Defaults to False.
    fixedKnots (array-like, optional):
        List of user-specified knots. Defaults to None.
    getAICc (bool, optional): 
        Return corrected Aikake information criterion? Defaults to False.
    natural (bool, optional): 
        Try a natural (no change in slope at the endpoints) cubic spline. Defaults to False. 

Returns:
    callable spline object
    AICc (optionally)
'''
```
```Python
getSplineErr(xdata, xspline, ydata, ydatae, num_knots=None, order=3, rand=False, fixedKnots=None, natural=False):
'''
Use getSpline to smooth mean and error bars. Create a spline-smooth band from that. 
'''
```
```Python
class TBSpline(xdata, ydata, edata=None, knots=None, order=3, naturalLike=False):
'''
A class that prepares a splrep and wraps it with splev.
'''
```
