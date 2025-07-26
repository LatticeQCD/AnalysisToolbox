latqcdtools.statistics.fitting
=============

```Python
do_fit(func, xdata, ydata, edata=None, start_params=None, priorval=None, priorsigma=None, algorithm='curve_fit', detailedInfo=False, show_results=False, **kwargs):
'''
Wrapper to fitter initialization and the fit in one step. See above for arguments. 
'''
```
```Python
save_func(func, filename, domain, args=(), func_err=None, args_err=(), grad=None, header=None, npoints=1000, **kwargs):
'''
'''
```
```Python
try_fit(func, xdata, ydata, edata=None, start_params=None, priorval=None, priorsigma=None, algorithms=['curve_fit', 'TNC', 'Powell', 'Nelder-Mead', 'nonlin'], detailedInfo=False, show_results=False, **kwargs):
'''
Wrapper to fitter initialization and the fit in one step. See above for arguments. For historical reasons
algorithms has no default values here. 
'''
```
```Python
unzipXYData(xydata):
'''Take a 2d xydata array, created by zipXYData, and extract xvalues and yvalues
for use inside of a function of two variables.

Args:
    xydata (np.ndarray): array of x,y coordinates [ (x1,y1), (x2,y1), ..., (x1,y2), ... ]

Returns:
    xvalues [x1, x2, ... , xN, x1, x2, ...],
    yavlues [y1, y1, ... , y1, y2, y2, ...]
'''
```
```Python
zipXYData(xdata, ydata):
'''Collect 1d xdata and ydata into an 2d xydata array. You can then use
unzipXYData inside of some func(xydata), which represents some f(x,y), to
separate the x part and y part.

Args:
    xdata (array-like)
    ydata (array-like)

Returns:
    np.ndarray: array of x,y coordinates [ (x1,y1), (x2,y1), ..., (x1,y2), ... ]
'''
```
```Python
class Fitter(func, xdata, ydata, edata=None, **kwargs):
'''
The :class:`Fitter`, contains all information necessary for fitting: The data, the function to be fitted, and
optional the data for the errors. There are different minimization algorithms available. Many of them need the
gradient or hessian of the chisquare. One way is to set the derivatives of the fitting function from outside.
The derivatives of the actual chisquare are then computed via error propagation. Another way is to use numerical
derivatives.

There are two ways to compute the derivatives of the chisqare numerically. Either compute the
numerical derivative of the whole chisquare (error_strat='hessian') or compute the derivatives of the fitting
function and use error propagation (error_strat='propagation'). The latter is the default case.
'''
```
