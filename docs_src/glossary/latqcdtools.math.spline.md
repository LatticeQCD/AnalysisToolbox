latqcdtools.math.spline
=============

```python
_averageWeights(b, X, Omega, lams, As, occam):
"""
Coefficients of the smoothing splines on the lam grid, which solve A_k c_k = X^T b, their weights,
proportional to their evidence, and their log evidences log p(y|lam), up to a constant that does not
depend on lam. b holds the whitened data and derivative constraints and X their rows, see
_smoothingSystem. The arrays As and occam have the lam grid as their first axis. The penalty is a Gaussian
prior on the coefficients with precision lam Omega, flat in the straight-line directions, which Omega does
not see. Integrating out c gives
    -2 log p(y|lam) = chi^2 + lam c^T Omega c + log det(A) - log det+(lam Omega) + const,
evaluated at the spline c, where det+ is the product of nonzero eigenvalues (REML). The last two terms are
-2 occam, see _smoothingAtLam.
"""
```
```python
_defined(vals, name, what):
"""
The entries of vals without nan. Warns how many of the what (e.g. samples) were dropped because the
observable name is undefined for them.
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
_fitRows(x, M, t, order, N, xd, od, ed):
"""
Rows X of the least-squares problem |X d - b|^2 of a spline of the given order with full knot vector t and
coefficients c = N d: the whitened design matrix M B N, followed by one row B^(od)(xd) N / ed for each
derivative constraint. The rows of b are the whitened data M y and the values dy/ed, see _Fitter._whitened.
"""
```
```python
_freePoint(idx, m):
"""
The m-th grid index, counting from 0, that is not in the configuration idx.
"""
```
```python
_kindSample(sample, fitter):
"""
Fit one bootstrap sample for bootSpline: ydata, followed by the values of the derivative constraints.
The bootstrap needs numbers, so pack what is needed to rebuild the spline into one flat array: the 
sample, the coefficients, and the number of parameters. bootSpline unpacks it. This lives at module 
level so that it can be pickled for nproc>1.
"""
```
```python
_knotChain(x, M, My, grid, lt, le, nu, kmax, natural, idx0, u, nequi, every, wmove, pbirth, pdeath, xd, od, ed):
"""
Metropolis-Hastings chain over knot configurations, starting from the configuration idx0, with u[it,:] the
6 uniform random numbers of update it. The first nequi updates are for thermalization; after them, every
every-th configuration is saved. xd, od, ed are the derivative constraints, see _knotFit. Returns a list
with the configuration, AIC, and coefficients of every saved configuration.
"""
```
```python
_knotFit(x, M, My, grid, idx, natural, xd, od, ed):
"""
Least-squares fit of the cubic regression spline with the configuration idx to the whitened data My, see
_lsqSystem. The last len(xd) entries of My are the whitened derivative constraints. Returns whether the fit
is determined, its AIC, and its B-spline coefficients.
"""
```
```python
_knotValid(idx, lt, le, nu):
"""
Is there a unique datum strictly inside each interval of the configuration idx?
"""
```
```python
_knotVector(xmin, xmax, interior, order=3):
"""
Full knot vector as in splrep: each endpoint repeated order+1 times around the interior knots.
"""
```
```python
_lamGrid(X, Omega, nu, nd, nlam):
"""
Grid of nlam values of lam, uniform in log lam, for averaging smoothing splines with nu unique x and nd
derivative constraints. tr(H) falls monotonically with lam and does not depend on y, so we bracket the
grid by the lam where tr(H) is close to its largest value (interpolation) and where it is close to 2
(straight line). The largest value is nu + nd, but at most the number of coefficients, nu + 2. This
makes the grid independent of the units of x and y. Each end is found by stepping in log lam by 1 (a
factor e in lam) from log lam = 0 and taking the last step at or beyond the target, so the grid is at most
a factor e wider at each end than the targets demand. The weights there are negligible, so the ends need
no more precision. If tr(H) does not reach a target before the normal equations become too
ill-conditioned to solve, which happens at very small or very large lam, the grid ends there instead.
Splines that far out carry practically no weight.
"""
```
```python
_lsqSystem(x, M, t, order, natural, xd, od, ed):
"""
Set up the least-squares problem |X d - b|^2 of a regression spline with full knot vector t, where
X = M B N, b = M y, and c = N d. For a natural spline, N imposes S''(a) = S''(b) = 0, see _naturalBasis;
otherwise N = 1. Each derivative constraint adds the row
B^(od)(xd) N / ed to X. The solution d solves the normal equations A d = X^T b with A = X^T X. Returns
N, X, A, and whether the fit cannot be solved: A is singular when the Schoenberg-Whitney conditions fail,
and too ill-conditioned when derivative constraints are much more precise than the data.
"""
```
```python
_modelSpreads(splines, w, xspl):
"""
Systematic errors of an average over models, e.g. knot configurations or values of lam: the weighted
spreads over the models, with weights w, of each observable of their splines, see _observables. Models
where an observable is undefined are dropped from its spread, see _spreadOf.
"""
```
```python
_naturalBasis(t):
"""
Matrix N with c = N d for the cubic splines with S''(a) = S''(b) = 0 and full knot vector t. The free
coefficients d are c_1,...,c_{nB-2}; the conditions fix c_0 and c_{nB-1}.
"""
```
```python
_nparams(k, natural, order=3)
```
```python
_observables(spline, xspl):
"""
The observables of a spline that bootSpline and posteriorSpline report, keyed by their names in the result:
the spline on xspl ('yspl'), the x on xspl where it is largest and its value there ('xmax', 'ymax'), and
its steepest inflection point ('xinflect', 'yinflect'), which is nan if there is none.
"""
```
```python
_parseDeriv(xdata, deriv, edata, maxorder):
"""
Interpret the soft derivative constraints deriv = (x, dy, de) or (x, dy, de, order): the order-th
derivative of the spline at x is dy with error de. Returns the arrays x, dy, de, order, which are empty
without deriv.
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
_smoothingAtLam(X, Omega, lam):
"""
The parts of the smoothing spline at lam that do not depend on y. Its coefficients solve the normal
equations A c = X^T b with A = X^T X + lam Omega, and A^-1 is also their posterior covariance, see
posteriorSpline.

Returns:
    A: the matrix of the normal equations
    Ainv: A^-1
    peff: effective number of parameters tr(H), where H = X A^-1 X^T maps b to the fit X c
    occam: the part of the log evidence that does not depend on y, see _averageWeights
"""
```
```python
_smoothingSystem(xdata, edata, xd, od, ed):
"""
Set up the least-squares problem of a cubic smoothing spline, which minimizes
(y - S)^T W (y - S) + sum_i ((S^(od_i)(xd_i) - dy_i)/ed_i)^2 + lam int S''^2, the middle term being the
derivative constraints. The minimizer over all twice-differentiable functions is a natural cubic spline
with a knot at every unique x, so we work in that B-spline basis. With S = B c and int S''^2 = c^T Omega c,
the first two terms are |X c - b|^2 with X = (M B, D) and b = (M y, dy/ed) stacked, where M is the
whitening matrix and D has the rows B^(od_i)(xd_i)/ed_i. None of this depends on y, dy, or lam.

Returns:
    t: full knot vector
    M: whitening matrix
    X: rows of the data and derivative constraints
    Omega: penalty matrix, int S''^2 = c^T Omega c
"""
```
```python
_spreadOf(vals, w):
"""
Weighted spread over models, with weights w, of their values vals of a quantity, which may be a number or
an array. Models where the quantity has a nan are dropped, and the remaining weights are renormalized.
Returns nan if it has a nan for every model.
"""
```
```python
_steepestInflection(spline, xspl):
"""
The inflection point (x,y) of spline on xspl where |spline'| is largest. Inflection points are the zeros
of the second derivative. Between knots, the second derivative is a polynomial of degree order-2, so its
zeros are computed exactly from its piecewise-polynomial form. Taking
the steepest one ignores the small wiggles a smoothed spline can have where the data are flat. A zero of
the second derivative at an end of xspl, as for a natural spline, does not count, so the first and last
grid intervals are skipped; numerically that zero is not exact. Pieces where the second derivative
vanishes identically give nan roots, which are skipped. Returns (nan,nan) if there is no inflection point.
"""
```
```python
_values(obs, name):
"""
The values of the observable name in the list obs of observables, see _observables.
"""
```
```python
bootSpline(xdata, ydata, edata, kind=Interpolating(), numb_samples=300, nsupport=301, seed=7271978, nproc=1, deriv=None) -> dict:
"""
Given xdata, ydata, edata, create a spline of the given kind. Use a Gaussian (parametric) bootstrap,
drawing ydata from normal distributions of width edata, to propagate uncertainties of the data into an
error band for the spline. Every bootstrap sample is fit with the same kind; e.g. for KnotAverage,
the chain is rerun in every sample. This bootstrap error is the statistical error. It does
not contain how much the models in an average (the knot configurations, or the values of lam)
disagree with each other, since the averaged spline can stay put while they disagree. For averaged kinds, the weighted
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
    deriv (tuple, optional):
        Soft derivative constraints, see getSpline. Their values are drawn from normal distributions of
        width de, independently of ydata, like the data. Defaults to None.

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
    For averaged kinds, KnotAverage and Smoothing(lam='average'), also:
        splinesModel, weightsModel: the spline of each model (knot configuration or value of lam) for
            the original data and its weight. Use modelSpread for the systematic error of your own
            observables.
"""
```
```python
getSpline(xdata, ydata, edata=None, kind=Interpolating(), deriv=None):
"""
Fit a spline of the given kind to the data.

Args:
    xdata (array-like)
    ydata (array-like)
    edata (array-like, optional):
        Error data, either a vector of errors or the covariance matrix of ydata. Defaults to None.
    kind (SplineKind, optional):
        Interpolating(), Regression(...), Smoothing(...), or KnotAverage(...). Defaults to Interpolating().
    deriv (tuple, optional):
        Soft constraints on derivatives of the spline, (x, dy, de) or (x, dy, de, order), each an array
        or a number: the order-th derivative of the spline at x is dy with error de. Each constraint adds
        ((S^(order)(x) - dy)/de)^2 to chi^2. order defaults to 1. Needs edata; not for Interpolating.
        Defaults to None.

Returns:
    TBSpline
"""
```
```python
getSplineErr(xdata, xspline, ydata, ydatae, kind=Interpolating()):
"""
Fit unweighted splines of the given kind to ydata-ydatae and ydata+ydatae, evaluate them at xspline,
and return their midpoint and half-difference. This gives a smooth band, but it is not a statistical
error propagation; for that use bootSpline.
"""
```
```python
modelSpread(res, f):
"""
Systematic error of a quantity f(spline) derived from a bootSpline result res: its weighted spread over
the models of the average (the knot configurations, or the values of lam) for the original data, the same
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
posteriorSpline(xdata, ydata, edata, kind=Smoothing(lam='average', nlam=100), numb_samples=300, nsupport=301, seed=7271978, deriv=None) -> dict:
"""
Smoothing spline and its Bayesian error band. The smoothing spline at fixed lam is the posterior mean
of a Bayesian model, in which the curvature penalty is a Gaussian prior on the B-spline coefficients.
Their posterior is Gaussian with covariance (B^T W B + D^T D + lam Omega)^-1, where D holds the rows of
the derivative constraints, see _smoothingSystem. With lam='average', we average over the lam grid,
weighted by the evidence, i.e. integrate lam out with a flat prior in log lam. This is not a bootstrap:
the data are used once, as measured, and the error band says how uncertain the true curve is, given the
data and the assumption that it is smooth. Unlike the bootstrap error of bootSpline, this includes the
uncertainty coming from the smoothing itself.

Args:
    xdata (array-like)
    ydata (array-like)
    edata (array-like):
        Errors of ydata, or their covariance matrix.
    kind (Smoothing, optional):
        Defaults to Smoothing(), i.e. the average over lam.
    numb_samples (int, optional):
        Number of samples drawn from the posterior. Defaults to 300.
    nsupport (int, optional):
        Number of points of the grid xspl, on which the spline is evaluated. Defaults to 301.
    seed (int, optional):
        Seed for the samples. Defaults to DEFAULTSEED.
    deriv (tuple, optional):
        Soft derivative constraints, see getSpline. They enter the likelihood like data, at their given
        values. Defaults to None.

Returns:
    dict: with entries
        xspl: grid
        yspl: posterior mean on xspl
        ysple: posterior standard deviation on xspl. ysple^2 = ysple_stat^2 + ysple_sys^2, where
            ysple_stat^2 = sum_k weights[k] Var[S|lam_k] is the average posterior variance at fixed
            lam and ysple_sys^2 the weighted variance of the posterior means over lam (law of total
            variance). These three are exact.
        splineMean: the posterior mean as a TBSpline
        splineSamples: numb_samples TBSplines drawn from the posterior. For a quantity derived from
            the spline, e.g. an integral or a derivative, take its spread over these samples. Do
            not add a systematic error on top; the spread already contains the uncertainty in lam.
        xmax, xmaxe: median and spread (dev_by_dist) over the samples of the x on xspl where the
            spline is largest. xmaxe_sys is the weighted spread over lam of the xmax of the
            posterior means, and xmaxe_stat = sqrt(xmaxe^2 - xmaxe_sys^2).
        ymax, ymaxe, xinflect, xinflecte, yinflect, yinflecte: the same for the spline at xmax and for the
            inflection point with the largest |slope|. Samples without an inflection point are dropped,
            and all are nan if no sample has one.
        lams, peffs, weights: the grid of lam, the effective numbers of parameters, and the weights.
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
class KnotAverage(kmax=None, ngrid=None, natural=False, nupdates=10000, nequi=1000, every=10, nchains=1, seed=7271978):
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
        Number of Metropolis-Hastings updates per chain after thermalization. Defaults to 10000.
    nequi (int, optional):
        Number of updates before those, discarded for thermalization. Defaults to 1000.
    every (int, optional):
        Save every this many updates after thermalization. Defaults to 10.
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
class Smoothing(lam='average', nlam=100):
"""
Cubic smoothing spline. It minimizes (y - S)^T W (y - S) + lam int S''^2, with W as for Regression.
For W = diag(w), this is the convention of scipy.interpolate.make_smoothing_spline. The minimizer is
a natural cubic spline with a knot at every unique x, so there are no knots to choose. With
lam='average', the result is instead an average of smoothing splines over nlam values of lam, uniform
in log lam, weighted by their evidence. The grid runs from (almost) interpolating the data to (almost)
a straight line. This needs edata. With derivative constraints, the result is the minimizer among the
cubic splines with a knot at every unique x. It is then no longer exactly natural: a constraint at an
end frees the curvature there, and one in between leaves a small curvature at the ends.

Args:
    lam (float or str, optional):
        Smoothing parameter > 0, or 'average'. Defaults to 'average'.
    nlam (int, optional):
        Number of values of lam for lam='average'. Defaults to 100.
"""
```
```python
class SplineKind():
"""
Base class of the kinds of spline. hasAIC is False for kinds without a meaningful AIC, i.e.
interpolation and averages over models. natural is True for kinds with zero curvature at the endpoints.
averaged is True for averages over models, whose fitters have details.
"""
```
```python
class TBSpline(xdata, ydata, tck, kind, edata=None, peff=None, deriv=None, lam=None, logEvidence=None):
"""
A fitted B-spline, as returned by getSpline, bootSpline, and posteriorSpline. It holds the data it was
fit to, the tck of the fit in the format of splrep, and the kind of spline, which also says whether
the spline is natural. Call it like a SciPy spline. For a spline averaged over knot configurations,
peff is the average number of parameters. deriv holds the derivative constraints of the fit as arrays
(x, dy, de, order), or None. For a smoothing spline, lam is the smoothing parameter, or 'average' for an
average over lam, and peff is the effective number of parameters. A smoothing spline at fixed lam also
knows its logEvidence.
"""
```
```python
class _Fitter(kind, xdata, edata, deriv=None):
"""
Base class of the fitters. Every fitter has fit(y, warn, dy), which returns a TBSpline; warn turns on
warnings about the fit, where a fitter has any, and dy replaces the values of the derivative
constraints, which bootSpline draws like data. Each subclass implements _fit(y, warn, dy).
"""
```
```python
class _InterpolatingFitter(kind, xdata, edata, deriv=None):
```
```python
class _KnotAverageFitter(kind, xdata, edata, deriv=None):
```
```python
class _RegressionFitter(kind, xdata, edata, deriv=None):
```
```python
class _SmoothingFitter(kind, xdata, edata, deriv=None):
```
