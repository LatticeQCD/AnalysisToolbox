# Curve fitting and splines


There are many ways to fit a curve. There are strategies that minimize $\chi^2/{\rm d.o.f.}$ 
for when you know the functional form ahead of time, splines for when you don't, and other methods. 
The AnalysisToolbox includes some routines that are helpful for this purpose.
By default the `Fitter` is [parallelized](../base/speedify.md) with `DEFAULTTHREADS`
processes. Set `nproc=1` if you want to turn off parallelization.
Splines are not parallelized in this way, but should be fast because they wrap SciPy methods.

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

The covariance matrix of the fit parameters are computed through error propagation of the covariance matrix
of the $y$-data. The error is then obtained by

$\sigma = \sqrt{\diag{\text{cov}}}$

In some codes, such as `gnuplot`, it is customary to multiply this error by a further factor $\chi^2/{\rm d.o.f.}$.
The intuition behind this is that the error will be increased if the fit is poor. This is okay if you would like to be
somewhat more conservative with your error bar, but it is strictly speaking not necessary. It also makes the error bar
more difficult to interpret clearly, i.e. if your input data and errors were well estimated, then it's not clear that
the true fit parameters will fall within one $\sigma$ of the estimators 67% of the time. The default behavior
is not to include this factor, but in case you would like it in your fits, for example because you are feeling
conservative, or for comparison with `gnuplot`, you can pass the option
`norm_err_chi2=True` to your `try_fit` or `do_fit` call.

## Splines

There are several methods in the toolbox to fit a 1D spline to some `xdata` and `ydata`.
These can be found in `latqcdtools.math.spline`. The basic method is `getSpline`
```Python
getSpline(xdata, ydata, num_knots=None, edata=None, order=3, rand=False, fixedKnots=None,
          getAICc=False, natural=False, seed=None)
```
By default this is a least-squares (regression) spline: a B-spline with fixed knots, fit to the data
by minimizing $\chi^2$ using `scipy.interpolate.splrep` with `task=-1`. (This is the same thing
`scipy.interpolate.LSQUnivariateSpline` does.) It does not in general pass through the data.
If you pass `edata`, each point is weighted by `1/edata`.
Here you specify the number of interior knots `num_knots` and the degree `order` of the spline polynomials.
Fewer knots give a stiffer curve. By default, the knots are evenly spaced in data index, so that each interval
between knots holds roughly the same number of data. You can specify `rand=True` to have it pick the knot
locations randomly; pass `seed` to make this reproducible. If you need to specify some knot locations
yourself, pass them as a list to `fixedKnots`. Note that `num_knots` includes these fixed knots in its counting;
hence if
```Python
len(fixedKnots)==num_knots
```
no knots will be generated automatically. With `getAICc=True`, which requires `edata`, `getSpline` also returns
the AICc of the fit, which you can use to compare different numbers of knots.

There also exists the option `natural=True` for natural splines, which have zero curvature at the endpoints.
This does two different things depending on whether you have errors:

- Without `edata`, you get an interpolating natural cubic spline from `scipy.interpolate.CubicSpline`.
  It has a knot at every data point and passes through all of them. Do not pass `num_knots` in this case.
- With `edata`, you get the least-squares spline described above, restricted to cubic splines with zero
  curvature at the endpoints. This removes two fit parameters, which reduces the variance of the spline
  near the ends of the data.

To propagate the errors of the data into an error band for the spline, use `bootSpline`, which
refits the spline to Gaussian bootstrap samples of the data:
```Python
res = bootSpline(xdata, ydata, edata, num_knots=num_knots)
plot_band(res['xspl'], res['yspl']-res['ysple'], res['yspl']+res['ysple'])
```
