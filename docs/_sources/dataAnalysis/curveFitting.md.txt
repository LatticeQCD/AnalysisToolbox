# Curve Fitting


There are many ways to fit a curve. There are strategies that minimize $\chi^2/{\rm d.o.f.}$ 
for when you know the functional form ahead of time, splines for when you don't, and other methods. 
The AnalysisToolbox includes some routines that are helpful for this purpose.

## $\chi^2$ minimization

In the module
```Python
import latqcdtools.statistics.fitting
``` 
one finds a `Fitter` class for carrying out fits. The `Fitter` class encapsulates all information
relevant to a fit, like its $x$-data, $y$-data, the fit form, and so on.
After constructing a `fitter` object, one can then use associated 
methods to try various kinds of fits. These are generally wrappers from `scipy.optimize`. 
An easy example is given in  `testing/fitting/simple_example.py`, shown below.

```Python
import numpy as np
import matplotlib.pyplot as plt
from latqcdtools.statistics.fitting import Fitter
from latqcdtools.base.logger import set_log_level

set_log_level('DEBUG')

print("\n Example of a simple 3-parameter quadratic fit.\n")

# Here we define our fit function. we pass it its independent variable followed by the fit parameters we are
# trying to determine.
def fit_func(x, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a*x**2 + b*x + c

xdata, ydata, edata = np.genfromtxt("wurf.dat", usecols=(0,2,3), unpack=True)

# We initialize our Fitter object. If expand = True, fit_func has to look like
#            func(x, a, b, *args)
#        otherwise it has to look like
#            func(x, params, *args).
fitter = Fitter(fit_func, xdata, ydata, expand = False)

# Here we try a fit, using the 'curve_fit' method, specifying the starting guesses for the fit parameters. Since
# ret_pcov = True, we will get back the covariance matrix as well.
res, res_err, chi_dof, pcov = fitter.try_fit(start_params = [1, 2, 3], algorithms = ['curve_fit'], ret_pcov = True)

print(" a , b,  c : ",res)
print(" ae, be, ce: ",res_err)
print("chi2/d.o.f.: ",chi_dof)
print("       pcov: \n",pcov,"\n")

fitter.plot_fit()
plt.show()
```

Supported fit algorithms include
- L-BFGS-B
- TNC
- Powell
- COBYLA
- SLSQP
which are essentially wrappers for `scipy` fit functions.
If one does not specify an algorithm, the `try_fit` method will attempt all of them and return the 
result of whichever one had the best $\chi^2/{\rm d.o.f.}$

**IMPORTANT: The functions that you pass to these fitting routines have to be able to handle arrays!** 
E.g. you pass it `[x0, x1, ..., xN]` and get back `[f(x0), f(x1), ..., f(xN)]`. It is written this 
way to force better performance; if it were a typical loop it would be slow. If you are having 
trouble figuring out how to write your function in a way to handle arrays, a good starting point 
can be to use [np.vectorize](https://numpy.org/doc/stable/reference/generated/numpy.vectorize.html).

## Splines

There are several methods in the toolbox to fit a 1D spline to some `xdata` and `ydata`.
These can be found in `latcqdtools.math.spline.py`. The basic method is `getSpline`
```Python
getSpline(xdata, ydata, num_knots, order=3, rand=False, fixedKnots=None)
```
This is essentially a wrapper for `scipy.interpolate.LSQUnivariateSpline`.
Here you specify how many knots `num_knots` you want and the order `order` of the spline.
By default, the spline will create a list of `num_knots` evenly spaced knots, but you can specify
`rand=True` to have it pick the knot locations randomly. If you need to specify some knot locations
yourself, this is accomplished by passing a list of specified knots to `fixedKnots`. Note that
`num_knots` includes these fixed knots in its counting; hence if
```Python
len(fixedKnots)==num_knots
```
no knots will be generated randomly.
