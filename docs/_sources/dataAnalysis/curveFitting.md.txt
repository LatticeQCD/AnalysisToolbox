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

The spline methods can be found in `latqcdtools.math.spline`. You choose what kind of spline you want by
passing a *kind* object. There are three kinds:

- `Interpolating()`: the natural cubic interpolating spline. It has a knot at every data point and
  passes through all of them.
- `Regression(nknots=None, knots=None, rand=False, seed=None, natural=False, order=3)`: a least-squares
  (regression) spline with fixed knots. It does not in general pass through the data.
- `KnotAverage(kmax=None, ngrid=None, natural=False, nupdates=10000, nburn=1000, every=10, nchains=1, seed=DEFAULTSEED)`:
  a cubic regression spline averaged over the number and positions of its knots.

The same kind can then be passed to either of the two main functions, which share their first four arguments:
```Python
spl = getSpline(xdata, ydata, edata=None, kind=None)          # one fit, returns a TBSpline
res = bootSpline(xdata, ydata, edata, kind=None, ...)         # bootstrap error band, returns a dict
```
As in the `Fitter`, a vector `edata` is read as errors and a matrix as the covariance matrix of `ydata`.
Every kind can take a covariance matrix. If you pass no kind, you get `Interpolating()`. The result of `getSpline` is a `TBSpline`, which you call like a
SciPy spline, e.g. `spl(x, der=1)` for its first derivative. `spl.get_nparams()` gives its number of parameters,
and `spl.get_AIC()`, which requires `edata`, its AIC, which you can use to compare e.g. different numbers of knots.
The AIC is not defined for `Interpolating`, which has no free parameters to compare, or for `KnotAverage`, which
already averages over models.

### Regression splines

`Regression` minimizes $\chi^2$ for B-splines with fixed knots. (Without correlations, this is the same thing
`scipy.interpolate.splrep` with `task=-1` and `scipy.interpolate.LSQUnivariateSpline` do.)
Here you specify the number of interior knots `nknots` and the degree `order` of the spline polynomials.
Fewer knots give a stiffer curve. By default, the knots are evenly spaced in data index, so that each interval
between knots holds roughly the same number of data. You can specify `rand=True` to have it pick the knot
locations randomly; pass `seed` to make this reproducible. Without `seed`, a seed is drawn once when the
kind is made, so all fits with this kind, e.g. in a bootstrap, use the same knots. If you need to specify some
knot locations yourself, pass them as a list to `knots`. Note that `nknots` includes these knots in its counting,
and defaults to `len(knots)`, in which case no knots are generated automatically.
With `natural=True` the spline has zero curvature at the endpoints. This removes two fit parameters, which
reduces the variance of the spline near the ends of the data.
```Python
spl = getSpline(xdata, ydata, edata, Regression(nknots=5, natural=True))
```

### Averaging over knots

Choosing the knots of a regression spline is a choice of model. To avoid making it, `KnotAverage` averages
regression splines over the number $k$ and the positions of their interior knots, so that knots can crowd
where the data demand structure. Each knot configuration gets the weight
$$
  p(k)\,p(t|k)\,e^{-{\rm AIC}/2},\qquad {\rm AIC}=\chi^2_{\rm min}+2p,
$$
with $p=k+4$ parameters ($k+2$ with `natural=True`). The knots lie on a grid of `ngrid` candidate positions,
by default fine enough to resolve every gap between data. $p(k)$ is uniform up to `kmax`, by default the largest
number of knots the data allow, and $p(t|k)$ is uniform over the configurations with $k$ knots on the grid.
Every interval between knots must contain at least one datum. The configurations are sampled by a Markov
chain (Metropolis–Hastings with birth, death, and move steps) of `nupdates` updates, of which the first
`nburn` are discarded and then every `every`-th is saved. The averaged spline is again a cubic spline, with a
knot at every grid point, and `get_nparams()` gives the average number of parameters. This needs `edata`.
```Python
spl = getSpline(xdata, ydata, edata, KnotAverage())
res = bootSpline(xdata, ydata, edata, KnotAverage(), nproc=8)
```
In `bootSpline` the chain is rerun for every bootstrap sample, and the spread over configurations for the
original data is added as a systematic error, see below; with the defaults, one chain takes about a second for 30 data, so use `nproc`. The chain is about
3 times faster if you call `numbaON()` before importing `latqcdtools.math.spline`, see [here](../base/speedify.md);
the result is the same either way. Each chain is seeded by
`seed`, the chain index, and the data, so results are reproducible and do not depend on `nproc`. The result
also contains diagnostics for the original data: the saved configurations `res['knotConfigs']` and their
counts `res['knotCounts']`, the histogram of $k$ in `res['kHistogram']`, and the acceptance rates. If the
histogram piles up at `kmax`, you get a warning. You also get a warning if a chain has not thermalized: the mean AIC
of the first 10% of its saved configurations is compared with that of the last 50% (Geweke test, errors from
batch means), and a difference of more than 3 standard errors, or fewer than 50 saved configurations, is flagged.
To check convergence further, compare runs with different `seed`, or longer chains.

### Error bands from the bootstrap

To propagate the errors of the data into an error band for the spline, use `bootSpline`. It draws
Gaussian bootstrap samples of the data with `bootstr_from_gauss` (correlated ones, if `edata` is a covariance
matrix) and fits each with the given kind. For a `Regression` the knots are the same for every sample.
Besides the band, it returns the fit to the original data in `res['splineMean']`, the bootstrap-level splines
in `res['splineBS']`, the location `res['xmax']` and height `res['ymax']` of the maximum and the steepest
inflection point `res['xinflect']`, `res['yinflect']`, all with errors, and, for kinds other than
`Interpolating` and averages, the AIC of the fit to the original data. Use `nproc` to fit the samples in parallel:
```Python
res = bootSpline(xdata, ydata, edata, Regression(nknots=5))
plot_band(res['xspl'], res['yspl']-res['ysple'], res['yspl']+res['ysple'])
```
The spread of the bootstrap samples is the statistical error. For an average over models, i.e. over
knot configurations, the weights are recomputed in every bootstrap sample, but this does not
capture how much the models disagree with each other: if the weights hardly change between samples, the
average hardly moves, however different the models are. For these kinds, the weighted spread over models for
the original data is added as a systematic error, as in a model average of fits. Every error comes with its
two parts, e.g. `ysple**2 = ysple_stat**2 + ysple_sys**2`, and likewise for `xmaxe`, `ymaxe`, `xinflecte`,
and `yinflecte`. For kinds that are not averages, the systematic part is 0. For a quantity of your own
derived from the spline, `res['splinesModel']` and `res['weightsModel']` hold the models for the original data
and their weights, and `modelSpread(res, f)` gives the weighted spread of `f(spline)` over them, to add in
quadrature to the spread of `f` over `res['splineBS']`. Being added in quadrature, the
total is somewhat conservative: the models also disagree partly because the noise affects them differently,
which the statistical error already contains.

Finally, `getSplineErr(xdata, xspline, ydata, ydatae, kind=None)` fits splines of the given kind (by default
`Interpolating()`) to `ydata-ydatae` and `ydata+ydatae` and returns their midpoint and half-difference on
`xspline`. This gives a smooth band, but it is not a statistical error propagation.
