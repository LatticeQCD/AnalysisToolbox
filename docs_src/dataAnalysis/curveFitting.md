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
          getAICc=False, natural=False, seed=None, lam=None)
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

If you pass a smoothing parameter `lam`, you get a cubic smoothing spline instead. It minimizes
$$
  \sum_i w_i\,\left(y_i-S(x_i)\right)^2 + \lambda\int dx\, S''(x)^2,
$$
with $w_i=1/\sigma_i^2$ from `edata` (or $w_i=1$ without `edata`), the same convention as
`scipy.interpolate.make_smoothing_spline`. The minimizer over all twice-differentiable functions is a natural
cubic spline with a knot at every unique $x$, so there are no knots to choose, and `num_knots`, `fixedKnots`,
and `rand` may not be passed. Unlike SciPy's version, several data may share one $x$.
Here $\lambda$ controls the smoothness: as $\lambda\to0$ the spline interpolates the data, and as
$\lambda\to\infty$ it becomes the weighted straight-line fit. Since the fit is linear in the data,
it has an effective number of parameters $p_{\rm eff}={\rm tr}\,H$, with $H$ the hat matrix mapping
the data to the fitted values. It runs from the number of unique $x$ down to 2, and it is what
`get_nparams()` returns and what `getAICc=True` uses. Note that $\lambda$ has units, since
$\int S''^2$ scales with the $x$ and $y$ ranges of your data, so a given value is only meaningful for one data set.
If your data are correlated, pass their covariance matrix as `edata`. As in the `Fitter`, a vector `edata`
is read as errors and a matrix as the covariance matrix. Then the first term above becomes
$(y-S)^T\,{\rm cov}^{-1}\,(y-S)$, and `getAICc=True` uses the covariance matrix as well. This is only
implemented for the smoothing spline.

The smoothing spline is the posterior mean of a Bayesian model, in which the penalty is a Gaussian prior
on the spline with precision $\lambda$ for its curvature, and flat for straight lines, which have no
curvature. The fit stores the log evidence $\log p(y|\lambda)$ of this model in its `logEvidence`
attribute, up to a constant that does not depend on $\lambda$. (This is the restricted likelihood, or
REML, of e.g. Wood, J. R. Stat. Soc. B 73, 3 (2011).)

To propagate the errors of the data into an error band for the spline, use `bootSpline`, which
refits the spline to Gaussian bootstrap samples of the data using `bootstr_from_gauss`. The knots are the
same for every sample. Besides the band, it returns the AICc of the fit to the original data, the
bootstrap-level splines in `res['splineBS']`, and the location of the maximum with its error. If your
data are correlated, pass their covariance matrix as `edata` to draw correlated samples. The fits themselves
are then weighted by the square root of its diagonal, unless you also pass `lam`: then each sample gets a
correlated smoothing spline. Use `nproc` to fit the samples in parallel:
```Python
res = bootSpline(xdata, ydata, edata, num_knots=num_knots)
plot_band(res['xspl'], res['yspl']-res['ysple'], res['yspl']+res['ysple'])
```

Choosing knots or $\lambda$ is a choice of model. To avoid making it, use `lam='average'`. This averages
smoothing splines over `nlam` values of $\lambda$, uniform in $\log\lambda$, weighted by their evidence,
i.e. it integrates $\lambda$ out with a flat prior in $\log\lambda$. Since the smoothing spline has a knot
at every unique $x$, there are no knots to choose either. The grid runs from $p_{\rm eff}$ close to the
number of unique $x$ (interpolation) to $p_{\rm eff}$ close to 2 (a straight line), so it does not depend
on the units of your data. In every bootstrap sample, the weights are recomputed, so the error `res['ysple']`
includes the uncertainty in $\lambda$. You also get the grid `res['lams']`, the corresponding `res['peffs']`,
`res['weights']`, and `res['splinesLam']`, the smoothing spline at each $\lambda$ for the original data.
`res['ysple_lam']` is the weighted spread of these splines. It is a diagnostic only: do not add it to
`res['ysple']`, which already contains the uncertainty in $\lambda$. If the ends of the grid carry noticeable
weight, you get a warning. At the large-$\lambda$ end this means the data are consistent with a straight line.
```Python
res = bootSpline(xdata, ydata, edata, lam='average')
plot_band(res['xspl'], res['yspl']-res['ysple'], res['yspl']+res['ysple'])
```

A bootstrap measures how much the spline would scatter if you repeated the experiment. It cannot see
the bias that the smoothing introduces, e.g. when it flattens a peak, and this matters most for derivatives.
`posteriorSpline` gives a Bayesian error band instead, which does include it. The smoothing spline at fixed
$\lambda$ is the posterior mean of a Bayesian model, in which the curvature penalty is a Gaussian prior
on the spline. `posteriorSpline` averages over $\lambda$ as above and draws `numb_samples` splines from
the resulting posterior: each sample picks a $\lambda$ with probability given by its weight, then a spline
from the Gaussian posterior at that $\lambda$. This is not a bootstrap; the data are used once, as measured.
`res['yspl']` is the posterior mean and `res['ysple']` the posterior standard deviation, which by the law
of total variance splits into `res['ysple_stat']` (the average posterior error at fixed $\lambda$) and
`res['ysple_sys']` (the spread over $\lambda$), with `ysple**2 = ysple_stat**2 + ysple_sys**2`.
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
