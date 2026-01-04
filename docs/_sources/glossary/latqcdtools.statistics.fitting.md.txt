latqcdtools.statistics.fitting
=============

```Python
do_fit(func, xdata, ydata, edata=None, start_params=None, priorval=None, priorsigma=None, algorithm='curve_fit', detailedInfo=False, show_results=False, **kwargs):
'''
Wrapper to fitter initialization and the fit in one step. See above for arguments. 
'''
```
```Python
save_func(func, filename, domain, args=(), func_err=None, args_err=(), grad=None, header=None, npoints=1000, **kwargs)
```
```Python
try_fit(func, xdata, ydata, edata=None, start_params=None, priorval=None, priorsigma=None, algorithms=['curve_fit', 'TNC', 'Powell', 'Nelder-Mead', 'nonlin'], detailedInfo=False, show_results=False, **kwargs):
'''
Wrapper to fitter initialization and the fit in one step. See above for arguments. For historical reasons
algorithms has no default values here. 
'''
```
```Python
unzipData(zipData):
'''Take a zipData array, created by zipData, and extract x-, y-, ... values 
for use inside of a function of two variables. Example in 2d:

Args:
    zipData (np.ndarray): In 2d, array of x,y coordinates [ (x1,y1), (x1,y2), ..., (x2,y1), ... ]

Returns: In 2d,
    xvalues [x1, x2, ..., xN, x1, x2, ..., xN, ...],
    yavlues [y1, y1, ..., y1, y2, y2, ..., y2, ...]
'''
```
```Python
zipData(*data):
'''Collect N sets of 1d data into an Nd zipData array. You can then use
unzipData inside of some func(zipData), which represents some f(x,y,...), to
separate the x-, y-, ... part.

Args:
    data: In 2d, would be e.g. xdata,ydata, both np.ndarrays

Returns:
    np.ndarray: In 2d, array of x,y coordinates [ (x1,y1), (x1,y2), ..., (x2,y1), ... ]
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
