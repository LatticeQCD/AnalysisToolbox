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
- `Smoothing(lam='average', nlam=100)`: a cubic smoothing spline with a curvature penalty, or an average of
  such splines over the smoothing parameter.
- `KnotAverage(kmax=None, ngrid=None, natural=False, nupdates=10000, nburn=1000, every=10, nchains=1, seed=DEFAULTSEED)`:
  a cubic regression spline averaged over the number and positions of its knots.

The same kind can then be passed to each of the three main functions, which share their first four arguments:
```Python
spl = getSpline(xdata, ydata, edata=None, kind=None)          # one fit, returns a TBSpline
res = bootSpline(xdata, ydata, edata, kind=None, ...)         # bootstrap error band, returns a dict
res = posteriorSpline(xdata, ydata, edata, kind=None, ...)    # Bayesian error band, Smoothing kinds only
```
As in the `Fitter`, a vector `edata` is read as errors and a matrix as the covariance matrix of `ydata`.
Every kind can take a covariance matrix. If you pass no kind, `getSpline` gives the interpolating spline
without `edata` and the average over $\lambda$ (see below) with `edata`, while `bootSpline` and
`posteriorSpline` use the average over $\lambda$. The result of `getSpline` is a `TBSpline`, which you call like a
SciPy spline, e.g. `spl(x, der=1)` for its first derivative. `spl.get_nparams()` gives its number of parameters,
and `spl.get_AICc()`, which requires `edata`, its AICc, which you can use to compare e.g. different numbers of knots.

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

### Smoothing splines

`Smoothing(lam=lam)` with a number `lam` gives a cubic smoothing spline. It minimizes
$$
  \sum_i w_i\,\left(y_i-S(x_i)\right)^2 + \lambda\int dx\, S''(x)^2,
$$
with $w_i=1/\sigma_i^2$ from `edata` (or $w_i=1$ without `edata`), the same convention as
`scipy.interpolate.make_smoothing_spline`. For a covariance matrix, the first term becomes
$(y-S)^T\,{\rm cov}^{-1}\,(y-S)$. The minimizer over all twice-differentiable functions is a natural
cubic spline with a knot at every unique $x$, so there are no knots to choose.
Unlike SciPy's version, several data may share one $x$.
Here $\lambda$ controls the smoothness: as $\lambda\to0$ the spline interpolates the data, and as
$\lambda\to\infty$ it becomes the weighted straight-line fit. Since the fit is linear in the data,
it has an effective number of parameters $p_{\rm eff}={\rm tr}\,H$, with $H$ the hat matrix mapping
the data to the fitted values. It runs from the number of unique $x$ down to 2, and it is what
`get_nparams()` returns and what `get_AICc()` uses. Note that $\lambda$ has units, since
$\int S''^2$ scales with the $x$ and $y$ ranges of your data, so a given value is only meaningful for one data set.

The smoothing spline is the posterior mean of a Bayesian model, in which the penalty is a Gaussian prior
on the spline with precision $\lambda$ for its curvature, and flat for straight lines, which have no
curvature. The fit stores the log evidence $\log p(y|\lambda)$ of this model in its `logEvidence`
attribute, up to a constant that does not depend on $\lambda$. (This is the restricted likelihood, or
REML, of e.g. Wood, J. R. Stat. Soc. B 73, 3 (2011).)

Choosing knots or $\lambda$ is a choice of model. To avoid making it, use `Smoothing()`, i.e.
`lam='average'`, which needs `edata`. This averages smoothing splines over `nlam` values of $\lambda$,
uniform in $\log\lambda$, weighted by their evidence, i.e. it integrates $\lambda$ out with a flat prior in
$\log\lambda$. The grid runs from $p_{\rm eff}$ close to the number of unique $x$ (interpolation) to
$p_{\rm eff}$ close to 2 (a straight line), so it does not depend on the units of your data. If the ends of the
grid carry noticeable weight, you get a warning. At the large-$\lambda$ end this means the data are consistent
with a straight line.

### Averaging over knots

The smoothing spline uses one $\lambda$ everywhere, so it cannot be stiff in one region and flexible in another.
This is a problem when the errors or the spacing of the data vary strongly. `KnotAverage` instead averages
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
In `bootSpline` the chain is rerun for every bootstrap sample, so the error includes the uncertainty in the
knots; with the defaults, one chain takes about a second for 30 data, so use `nproc`. The chain is about
3 times faster if you call `numbaON()` before importing `latqcdtools.math.spline`, see [here](../base/speedify.md);
the result is the same either way. Each chain is seeded by
`seed`, the chain index, and the data, so results are reproducible and do not depend on `nproc`. The result
also contains diagnostics for the original data: the saved configurations `res['knotConfigs']` and their
counts `res['knotCounts']`, the histogram of $k$ in `res['kHistogram']`, and the acceptance rates. If the
histogram piles up at `kmax`, you get a warning. `res['ysple_knots']`, the spread over configurations, is
a diagnostic only, like `res['ysple_lam']`. You also get a warning if a chain has not thermalized: the mean AIC
of the first 10% of its saved configurations is compared with that of the last 50% (Geweke test, errors from
batch means), and a difference of more than 3 standard errors, or fewer than 50 saved configurations, is flagged.
To check convergence further, compare runs with different `seed`, or longer chains.

### Error bands from the bootstrap

To propagate the errors of the data into an error band for the spline, use `bootSpline`. It draws
Gaussian bootstrap samples of the data with `bootstr_from_gauss` (correlated ones, if `edata` is a covariance
matrix) and fits each with the given kind. For a `Regression` the knots are the same for every sample.
Besides the band, it returns the fit to the original data in `res['splineMean']`, the bootstrap-level splines
in `res['splineBS']`, the location of the maximum with its error, and, for kinds other than `Interpolating` and
averages, the AICc of the fit to the original data. Use `nproc` to fit the samples in parallel:
```Python
res = bootSpline(xdata, ydata, edata, Regression(nknots=5))
plot_band(res['xspl'], res['yspl']-res['ysple'], res['yspl']+res['ysple'])
```
For the average over $\lambda$, the weights are recomputed in every bootstrap sample, so the error `res['ysple']`
includes the uncertainty in $\lambda$. You also get the grid `res['lams']`, the corresponding `res['peffs']`,
`res['weights']`, and `res['splinesLam']`, the smoothing spline at each $\lambda$ for the original data.
`res['ysple_lam']` is the weighted spread of these splines. It is a diagnostic only: do not add it to
`res['ysple']`, which already contains the uncertainty in $\lambda$.
```Python
res = bootSpline(xdata, ydata, edata, Smoothing())
plot_band(res['xspl'], res['yspl']-res['ysple'], res['yspl']+res['ysple'])
```

### Error bands from the posterior

A bootstrap measures how much the spline would scatter if you repeated the experiment. It cannot see
the bias that the smoothing introduces, e.g. when it flattens a peak, and this matters most for derivatives.
`posteriorSpline` gives a Bayesian error band instead, which does include it. It takes `Smoothing` kinds only.
The smoothing spline at fixed $\lambda$ is the posterior mean of a Bayesian model, in which the curvature
penalty is a Gaussian prior on the spline. For the average over $\lambda$, `posteriorSpline` draws
`numb_samples` splines from the resulting posterior: each sample picks a $\lambda$ with probability given by
its weight, then a spline from the Gaussian posterior at that $\lambda$. This is not a bootstrap; the data are
used once, as measured. `res['yspl']` is the posterior mean and `res['ysple']` the posterior standard deviation,
which by the law of total variance splits into `res['ysple_stat']` (the average posterior error at fixed
$\lambda$) and `res['ysple_sys']` (the spread over $\lambda$), with `ysple**2 = ysple_stat**2 + ysple_sys**2`.
For a quantity derived from the spline, e.g. an integral or a derivative, compute it for each of the samples
in `res['splineSamples']` and take their spread. Do not add a systematic error on top of that, since the
samples already contain the uncertainty in $\lambda$. In a test with a smooth function, both error bands
covered the truth about as often as they should, but for its derivative the posterior band was conservative
(89% instead of 68%), while the bootstrap band was close to nominal.
```Python
res = posteriorSpline(xdata, ydata, edata)
plot_band(res['xspl'], res['yspl']-res['ysple'], res['yspl']+res['ysple'])
derivs = []
for spl in res['splineSamples']:
    derivs.append(spl(res['xspl'], der=1))
dm, de = std_median(derivs), dev_by_dist(derivs)
```
Keep in mind that this is still a model: it assumes the same smoothness everywhere, and the smoothing
spline has zero curvature at the ends of the data, which biases it there if your data curve strongly.

Finally, `getSplineErr(xdata, xspline, ydata, ydatae, kind=None)` fits splines of the given kind (by default
`Interpolating()`) to `ydata-ydatae` and `ydata+ydatae` and returns their midpoint and half-difference on
`xspline`. This gives a smooth band, but it is not a statistical error propagation.
