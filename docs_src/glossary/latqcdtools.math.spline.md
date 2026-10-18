latqcdtools.math.spline
=============

```python
_batchMeanError(a, nbatch=5):
"""
Error of the mean of the time series a from the spread of the means of nbatch consecutive batches.
"""
```
```python
_checkEdata(xdata, edata)
```
```python
_dataSeed(seed, ichain, y):
"""
Seed for one chain from the seed of the kind, the chain index, and the bits of the data.
"""
```
```python
_even_knots(xdata, nknots):
"""
Return a list of nknots knots, evenly spaced in data index, so that each interval holds
roughly the same number of data. Each knot is the midpoint of two neighboring unique data, or a datum
itself when the index lands exactly on one.
"""
```
```python
_kindSample(sample, fitter):
"""
Fit one bootstrap sample of shape (1,ndata) for bootSpline. The bootstrap needs numbers, so pack what
is needed to rebuild the spline into one flat array: yBS, the coefficients, and the number of parameters.
_unpackSample undoes this. This lives at module level so that it can be pickled for nproc>1.
"""
```
```python
_knotVector(xmin, xmax, interior, order=3):
"""
Full knot vector as in splrep: each endpoint repeated order+1 times around the interior knots.
"""
```
```python
_lsqSystem(x, M, t, order, natural):
"""
Set up the least-squares problem |M B N d - M y|^2 of a regression spline with full knot vector t by QR, 
where c = N d. For a natural spline, the conditions S''(a) = S''(b) = 0 are C c = 0 and N is a basis of 
the null space of C; otherwise N = 1. Returns N, Q, R, and whether the fit is underdetermined, i.e. 
whether the Schoenberg-Whitney conditions fail.
"""
```
```python
_maxOnGrid(ys, xspl):
"""
For curves ys on xspl, one per row, the x on xspl where each is largest, and its value there.
"""
```
```python
_medianBiasCorrect(f0, fBS):
"""
Bias-corrected central value 2 f0 - median(fBS), with f0 the observable for the original data and fBS
its bootstrap samples along axis 0. This is _biasCorrect of the bootstrap with the median in place of
the mean.
"""
```
```python
_modelSpreads(splines, ys, w, xspl):
"""
Systematic errors of an average over models (knot configurations): the weighted spreads over the
models, with weights w, of each model's spline on xspl (ys holds them), its maximum (x and y on xspl), and
its steepest inflection point. Models without an inflection point are dropped from that spread, see
_spreadOf.
"""
```
```python
_nparams(k, natural):
"""
Number of parameters of a cubic regression spline with k interior knots. Natural boundary conditions
remove 2.
"""
```
```python
_parseEdata(edata):
"""
Interpret edata like the Fitter does: a vector holds errors and a matrix is the covariance matrix.
Returns the weights 1/sigma and the covariance matrix, which is None unless edata is a matrix.
"""
```
```python
_random_knots(xdata, nknots, SEED=None):
"""
Return a list of nknots randomly placed knots. Draws a random subset of nknots+1 data and applies
_even_knots to it.
"""
```
```python
_setErrors(res, name, ename, f0, fBS, sys):
"""
Put the bias-corrected central value of an observable into res[name], and its statistical error over
the bootstrap samples fBS, its systematic error sys, and their sum in quadrature into res[ename+'_stat'],
res[ename+'_sys'], and res[ename].
"""
```
```python
_spreadOf(splines, w, f):
"""
Weighted spread over the models splines, with weights w, of f(spline), which may be a number or an array.
Models where f has a nan are dropped, and the remaining weights are renormalized. Returns nan if f has a
nan for every model.
"""
```
```python
_steepestInflection(spline, xspl):
"""
The inflection point (x,y) of spline on xspl where |spline'| is largest. Inflection points are where the
second derivative changes sign; each root is found between the two grid points that bracket it. Taking
the steepest one ignores the small wiggles a smoothed spline can have where the data are flat. A zero of
the second derivative at an end of xspl, as for a natural spline, does not count, so the first and last
grid intervals are skipped; numerically that zero is not exact. Returns (nan,nan) if there is no
inflection point.
"""
```
```python
_unpackSample(sample, ndata):
"""
yBS, the coefficients, and the number of parameters from an array packed by _kindSample.
"""
```
```python
_weightedSpread(w, a):
"""
Spread sqrt(sum_k w_k (a_k - abar)^2) of a along axis 0, with abar = sum_k w_k a_k and normalized w.
"""
```
```python
_whitening(edata, ndata):
"""
Whitening matrix M with M^T M = W, so that a chi^2 (y - S)^T W (y - S) becomes |M y - M S|^2. For
errors, M = diag(1/edata). For a covariance matrix, M is the whiteningMatrix of the toolbox. Without
edata, M = 1.
"""
```
```python
bootSpline(xdata, ydata, edata, kind=None, numb_samples=300, nsupport=301, seed=7271978, nproc=1) -> dict:
"""
Given xdata, ydata, edata, create a spline of the given kind. Use a Gaussian (parametric) bootstrap,
drawing ydata from normal distributions of width edata, to propagate uncertainties of the data into an
error band for the spline. Every bootstrap sample is fit with the same kind; e.g. for KnotAverage,
the chain is rerun in every sample. This bootstrap error is the statistical error. It does
not contain how much the models in an average (the knot configurations) disagree with each
other, since the averaged spline can stay put while they disagree. For averaged kinds, the weighted
spread over models for the original data is therefore the systematic error, and every error below is
the statistical and systematic errors added in quadrature. If edata is the covariance matrix of ydata,
the samples are drawn from a multivariate normal, and the fits use the covariance matrix as well.

Args:
    xdata (array-like)
    ydata (array-like)
    edata (array-like):
        Errors of ydata, or their covariance matrix.
    kind (SplineKind, optional):
        Kind of spline. Defaults to Interpolating().
    numb_samples (int, optional):
        Number of bootstrap samples. Defaults to 300.
    nsupport (int, optional):
        Number of points of the grid xspl, on which the spline is evaluated. Defaults to 301.
    seed (int, optional):
        Seed for the samples. Defaults to DEFAULTSEED.
    nproc (int, optional):
        Number of processes, passed to bootstr_from_gauss. Defaults to 1.

Returns:
    dict: with entries
        xspl: grid
        yspl, ysple: central value and error on xspl. The central value is 2 f - median(f*), with f from
            the original data and f* from the samples. ysple^2 = ysple_stat^2 + ysple_sys^2, where
            ysple_stat is dev_by_dist of the samples and ysple_sys the weighted spread over models for
            the original data, which is 0 unless the kind is averaged.
        xmax, xmaxe, ymax, ymaxe: the same for the x on xspl where the spline is largest and the spline
            there, each error with its _stat and _sys parts
        xinflect, xinflecte, yinflect, yinflecte: the same for the inflection point with the largest |slope|.
            Samples and models without an inflection point are dropped, and all are nan if splineMean
            has none.
        splineMean: the spline fit to the original data, not an average over samples
        splineBS: the spline of each bootstrap sample
        AIC: AIC of splineMean, unless the kind is Interpolating or averaged
    For KnotAverage also, for the original data:
        knotConfigs, knotCounts: the distinct saved knot configurations and how often each was saved.
        kHistogram: how often each number of knots 0,...,kmax was saved.
        acceptance: acceptance rates of birth, death, and move proposals.
        ngrid, kmax: the number of candidate knots and the largest number of knots.
    For averaged kinds also:
        splinesModel, weightsModel: the spline of each model (knot configuration) for the original
            data and its weight. Use modelSpread for the systematic error of your own observables.
"""
```
```python
getSpline(xdata, ydata, edata=None, kind=None):
"""
Fit a spline of the given kind to the data.

Args:
    xdata (array-like)
    ydata (array-like)
    edata (array-like, optional):
        Error data, either a vector of errors or the covariance matrix of ydata. Defaults to None.
    kind (SplineKind, optional):
        Interpolating(), Regression(...), or KnotAverage(...). Defaults to Interpolating().

Returns:
    TBSpline
"""
```
```python
getSplineErr(xdata, xspline, ydata, ydatae, kind=None):
"""
Fit unweighted splines of the given kind to ydata-ydatae and ydata+ydatae, evaluate them at xspline,
and return their midpoint and half-difference. This gives a smooth band, but it is not a statistical
error propagation; for that use bootSpline. kind defaults to Interpolating().
"""
```
```python
modelSpread(res, f):
"""
Systematic error of a quantity f(spline) derived from a bootSpline result res: its weighted spread over
the models of the average (the knot configurations) for the original data, the same
systematic error bootSpline gives its own observables. Add it in quadrature to the statistical error,
i.e. the spread of f over res['splineBS']. Models where f is nan are dropped. For kinds that are not
averages, there are no models, and the result is 0.

Args:
    res (dict): result of bootSpline
    f (func): takes a TBSpline and returns a number or an array

Returns:
    float or np.ndarray: systematic error of f
"""
```
```python
class Interpolating():
"""
Natural cubic interpolating spline, with a knot at every datum. It passes through the data, so edata
are not used in the fit.
"""
```
```python
class KnotAverage(kmax=None, ngrid=None, natural=False, nupdates=10000, nburn=1000, every=10, nchains=1, seed=7271978):
"""
Cubic regression spline averaged over the number and positions of its interior knots, with AIC weights
p(k) p(t|k) exp(-AIC/2), AIC = chi^2_min + 2 p. The knots lie on a grid of ngrid candidate positions 
strictly inside the data range; p(k) is uniform on 0,...,kmax and p(t|k) uniform over the configurations 
with k knots on the grid. Configurations are sampled by Metropolis-Hastings with birth, death, and move 
steps, and every configuration needs at least one unique datum strictly between neighboring knots and 
between the ends of the data and the outermost knots. The average is again a cubic spline, with a knot 
at every grid point. This needs edata. In bootSpline the chain is rerun for every bootstrap sample.

Args:
    kmax (int, optional):
        Largest number of interior knots. Defaults to the largest number the data allow.
    ngrid (int, optional):
        Number of candidate knot positions. Defaults to max(100, 4 (xmax-xmin)/(smallest gap between 
        unique x)), so that the grid resolves every gap between data.
    natural (bool, optional):
        Zero curvature at the endpoints, which removes 2 parameters. Defaults to False.
    nupdates (int, optional):
        Number of Metropolis-Hastings updates per chain, including nburn. Defaults to 10000.
    nburn (int, optional):
        Number of updates discarded for thermalization. Defaults to 1000.
    every (int, optional):
        Save every this many updates after nburn. Defaults to 10.
    nchains (int, optional):
        Number of independent chains, whose saved configurations are pooled. Defaults to 1.
    seed (int, optional):
        Seed of the chains. Each chain is seeded by seed, its index, and the data, so that every bootstrap 
        sample gets its own reproducible chain. Defaults to DEFAULTSEED.
"""
```
```python
class Regression(nknots=None, knots=None, rand=False, seed=None, natural=False, order=3):
"""
Least-squares (regression) spline with fixed interior knots. It minimizes (y - S)^T W (y - S), with
W = diag(1/edata^2), or edata^-1 if edata is a covariance matrix, or 1 without edata. By default the
knots are evenly spaced in data index, so that each interval between knots holds roughly the same
number of data.

Args:
    nknots (int, optional):
        Number of interior knots, including knots. Defaults to len(knots).
    knots (list, optional):
        Interior knots you choose yourself. These count toward nknots. Defaults to None.
    rand (bool, optional):
        Place the remaining knots randomly? Defaults to False.
    seed (int, optional):
        Seed for rand. Defaults to None, which draws a seed once, so that every fit with this kind
        uses the same knots.
    natural (bool, optional):
        Zero curvature at the endpoints, which removes 2 parameters. Defaults to False.
    order (int, optional):
        Degree of the spline polynomials (SciPy's k). Defaults to 3.
"""
```
```python
class SplineKind():
"""
Base class of the kinds of spline. hasAIC is False for kinds without a meaningful AIC, i.e.
interpolation and averages over models. natural is True for kinds with zero curvature at the endpoints.
"""
```
```python
class TBSpline(xdata, ydata, tck, kind, edata=None, peff=None):
"""
A fitted B-spline, as returned by getSpline and bootSpline. It holds the data it was
fit to, the tck of the fit in the format of splrep, and the kind of spline, which also says whether
the spline is natural. Call it like a SciPy spline. For a spline averaged over knot configurations,
peff is the average number of parameters.
"""
```
```python
class _Fitter(kind, xdata, edata):
"""
Base class of the fitters. Every fitter has fit(y, warn), which returns a TBSpline; warn turns on
warnings about the fit, where a fitter has any. Fitters of averaged kinds also have details(y, xspl),
which returns the models of the average for the data y and their weighted spreads.
"""
```
```python
class _InterpolatingFitter(kind, xdata, edata):
```
```python
class _KnotAverageFitter(kind, xdata, edata):
```
```python
class _RegressionFitter(kind, xdata, edata):
```
