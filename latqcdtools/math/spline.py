#
# spline.py
#
# D. Clarke, Claude Code
#
# Splines with control over knot placement, endpoints, and smoothing, and their error bands. What kind of
# spline you get is set by a kind object: Interpolating, Regression, Smoothing, or KnotAverage. Pass it to
# getSpline for a single fit, to bootSpline for a bootstrap error band, or, for Smoothing, to
# posteriorSpline for a Bayesian one.
#

import numpy as np
from scipy.interpolate import BSpline, make_interp_spline, splev
from scipy.linalg import block_diag
from scipy.optimize import brentq
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import AIC, whiteningMatrix, std_median, dev_by_dist, getModelWeights, \
    parseEdata, checkEdata, weightedSpread
from latqcdtools.base.check import checkType, checkEqualLengths
from latqcdtools.base.initialize import TBRNG, DEFAULTSEED
from latqcdtools.math.math import isMatrix, quadrature, logDet
from latqcdtools.statistics.bootstr import bootstr_from_gauss, biasCorrect
from latqcdtools.statistics.autocorrelation import gewekeZ


def _even_knots(xdata, nknots):
    """
    Return a list of nknots knots, evenly spaced in data index, so that each interval holds
    roughly the same number of data. Each knot is the midpoint of two neighboring unique data, or a datum
    itself when the index lands exactly on one.
    """
    if len(xdata)<nknots:
        logger.TBRaise('number of data < number of knots')
    # sorted, without duplicate x values
    flat_xdata = np.unique(np.asarray(xdata))
    jump_step = (len(flat_xdata) - 1) / (nknots + 1)
    knots = []
    for i in range(1, nknots + 1):
        x_lower = flat_xdata[int(np.floor(i * jump_step))]
        x_upper = flat_xdata[int(np.ceil(i * jump_step))]
        knots.append((x_lower + x_upper) / 2)
    return knots


def _random_knots(xdata, nknots, SEED=None):
    """
    Return a list of nknots randomly placed knots. Draws a random subset of nknots+1 data and applies
    _even_knots to it.
    """
    rng = TBRNG(SEED)
    flat_xdata = np.sort(np.asarray(xdata))
    nsample = nknots+1
    sample_xdata = rng.choice(flat_xdata,nsample,replace=False)
    # Retry if too many data points are removed by np.unique
    while len(np.unique(sample_xdata)) < nsample:
        sample_xdata = rng.choice(flat_xdata,nsample,replace=False)
    return _even_knots(sample_xdata, nknots)


def _parseDeriv(xdata, deriv, edata, maxorder):
    """
    Interpret the soft derivative constraints deriv = (x, dy, de) or (x, dy, de, order): the order-th
    derivative of the spline at x is dy with error de. Returns the arrays x, dy, de, order, which are empty
    without deriv.
    """
    if deriv is None:
        return np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0,dtype=np.int64)
    if len(deriv) == 3:
        xd, dy, de = deriv
        order = np.ones(np.size(xd),dtype=np.int64)
    elif len(deriv) == 4:
        xd, dy, de, order = deriv
    else:
        logger.TBRaise('deriv must be (x, dy, de) or (x, dy, de, order).')
    xd    = np.atleast_1d(np.asarray(xd,dtype=float))
    dy    = np.atleast_1d(np.asarray(dy,dtype=float))
    de    = np.atleast_1d(np.asarray(de,dtype=float))
    order = np.atleast_1d(np.asarray(order,dtype=np.int64))
    checkEqualLengths(xd,dy,de,order)
    if edata is None:
        logger.TBRaise('deriv needs edata, since the derivative constraints are weighted against the data.')
    if np.any(de <= 0):
        logger.TBRaise('Need deriv errors de > 0. de =',de)
    if np.any(order < 1) or np.any(order > maxorder):
        logger.TBRaise('Need 1 <= deriv order <=',maxorder,'. order =',order)
    if np.any(xd < np.min(xdata)) or np.any(xd > np.max(xdata)):
        logger.TBRaise('deriv positions must lie inside the range of xdata. x =',xd)
    return xd, dy, de, order


# Largest condition number of the matrix A of the normal equations A c = r that we accept. Solving them loses up to
# log10 cond(A) of the 16 digits of double precision.
_CONDMAX = 1e14


def _fitRows(x, M, t, order, N, xd, od, ed):
    """
    Rows X of the least-squares problem |X d - b|^2 of a spline of the given order with full knot vector t and
    coefficients c = N d: the whitened design matrix M B N, followed by one row B^(od)(xd) N / ed for each
    derivative constraint. The rows of b are the whitened data M y and the values dy/ed, see _Fitter._whitened.
    """
    nB   = len(t) - order - 1
    B    = BSpline.design_matrix(x, t, order).toarray()
    rows = [M @ B @ N]
    for i in range(len(xd)):
        dB = BSpline(t, np.eye(nB), order).derivative(od[i])(xd[i])
        rows.append(np.atleast_2d(dB @ N)/ed[i])
    return np.vstack(rows)


def _smoothingSystem(xdata, edata, xd, od, ed):
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
    x = np.asarray(xdata,dtype=float)
    u = np.unique(x)
    t = _knotVector(u[0], u[-1], u[1:-1])
    nB = len(t) - 4
    M  = whiteningMatrix(parseEdata(edata, len(x)))
    # B_j'' is linear between neighboring knots, so B_i'' B_j'' is quadratic there, and 2-point
    # Gauss-Legendre quadrature on each interval gives int S''^2 = sum_q wq S''(xq)^2 exactly.
    half  = np.diff(u)/2
    mid   = (u[:-1]+u[1:])/2
    xq    = np.r_[mid - half/np.sqrt(3), mid + half/np.sqrt(3)]
    wq    = np.r_[half, half]
    L     = np.sqrt(wq)[:,None]*BSpline(t, np.eye(nB), 3).derivative(2)(xq)
    Omega = L.T @ L
    X     = _fitRows(x, M, t, 3, np.eye(nB), xd, od, ed)
    return t, M, X, Omega


def _smoothingAtLam(X, Omega, lam):
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
    A = X.T @ X + lam*Omega
    if np.linalg.cond(A) > _CONDMAX:
        logger.TBRaise('The smoothing spline at lam =',lam,'cannot be solved accurately, since its normal equations',
                       'have condition number',np.linalg.cond(A),'. This happens for very large lam, and when derivative',
                       'constraints are much more precise than the data; then give them larger errors.')
    Ainv = np.linalg.inv(A)
    # tr(X A^-1 X^T) = tr(A^-1 X^T X). Solving is more accurate than multiplying by A^-1.
    peff = np.trace(np.linalg.solve(A, X.T @ X))
    # Omega misses the straight lines, so its rank is nB-2.
    nB    = len(A)
    occam = -0.5*logDet(A) + 0.5*(nB-2)*np.log(lam)
    return A, Ainv, peff, occam


def _lamGrid(X, Omega, nu, nd, nlam):
    """
    Grid of nlam values of lam, uniform in log lam, for averaging smoothing splines with nu unique x and nd
    derivative constraints. tr(H) falls monotonically with lam and does not depend on y, so we bracket the
    grid by the lam where tr(H) is close to its largest value (interpolation) and where it is close to 2
    (straight line). The largest value is nu + nd, but at most the number of coefficients, nu + 2. This
    makes the grid independent of the units of x and y. If tr(H) does not reach a target before the normal
    equations become too ill-conditioned to solve, which happens at very small or very large lam, the grid
    ends there instead. Splines that far out carry practically no weight.
    """
    if nu < 3:
        logger.TBRaise('Averaging over lam needs at least 3 unique x.')
    peffmax = min(nu + nd, nu + 2)

    def peffDiff(loglam,target):
        _, _, peff, _ = _smoothingAtLam(X, Omega, np.exp(loglam))
        return peff - target

    def solvable(loglam):
        return np.linalg.cond(X.T @ X + np.exp(loglam)*Omega) <= _CONDMAX

    # Bracket each end in steps of a factor e in lam.
    loglams = []
    for target in [peffmax-0.1, 2.01]:
        lo, hi = -1., 1.
        while peffDiff(lo,target) < 0 and solvable(lo-1):
            lo -= 1
        while peffDiff(hi,target) > 0 and solvable(hi+1):
            hi += 1
        if peffDiff(lo,target) < 0:
            loglams.append(lo)
        elif peffDiff(hi,target) > 0:
            loglams.append(hi)
        else:
            loglams.append(brentq(peffDiff,lo,hi,args=(target,)))
    return np.exp(np.linspace(loglams[0],loglams[1],nlam))


def _averageWeights(b, X, Omega, lams, As, occam):
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
    Xb = X.T @ b
    c, logEv = [], []
    for k in range(len(lams)):
        ck = np.linalg.solve(As[k], Xb)
        r  = b - X @ ck
        c.append(ck)
        logEv.append(-0.5*(r @ r + lams[k] * ck @ Omega @ ck) + occam[k])
    logEv = np.array(logEv)
    # The log evidence plays the role of -IC/2.
    return getModelWeights(-2*logEv), logEv, np.array(c)


def _nparams(k, natural, order=3):
    nparams = k + order + 1
    if natural:
        nparams -= 2
    return nparams


class TBSpline:

    """
    A fitted B-spline, as returned by getSpline, bootSpline, and posteriorSpline. It holds the data it was
    fit to, the tck of the fit in the format of splrep, and the kind of spline, which also says whether
    the spline is natural. Call it like a SciPy spline. For a spline averaged over knot configurations,
    peff is the average number of parameters. deriv holds the derivative constraints of the fit as arrays
    (x, dy, de, order), or None. For a smoothing spline, lam is the smoothing parameter, or 'average' for an
    average over lam, and peff is the effective number of parameters. A smoothing spline at fixed lam also
    knows its logEvidence.
    """

    def __init__(self,xdata,ydata,tck,kind,edata=None,peff=None,deriv=None,lam=None,logEvidence=None):
        self.xspline     = np.copy(xdata)
        self.yspline     = np.copy(ydata)
        # Store tck in the format of splrep, which pads the coefficients with order+1 zeros. SciPy's splder and
        # splantider need that padding.
        t, c, k = tck
        c = np.asarray(c,dtype=float)
        if len(c) == len(t) - k - 1:
            c = np.r_[c, np.zeros(k+1)]
        self.tck         = (t, c, k)
        self.kind        = kind
        self.edata       = edata
        self.peff        = peff
        self.deriv       = deriv
        self.lam         = lam
        self.logEvidence = logEvidence
        # The covariance matrix of the data
        self.cov = None
        if edata is not None:
            self.cov = parseEdata(edata)

    @property
    def natural(self):
        return self.kind.natural

    def get_nparams(self):
        """
        Number of free parameters. Natural boundary conditions remove 2. For a smoothing spline this
        is the effective number tr(H), which is not an integer. For an average over lam or over knot
        configurations, it is the weighted average.
        """
        if self.peff is not None:
            return self.peff
        # The full knot vector repeats each end order+1 times.
        k = len(self.get_knots()) - 2*(self.get_order() + 1)
        return _nparams(k, self.natural, self.get_order())

    def get_AIC(self):
        """
        AIC of the fit, using the covariance matrix if edata is one. Requires edata, and a kind with a
        meaningful AIC. Derivative constraints count like data in -2 log L.
        """
        if not self.kind.hasAIC:
            logger.TBRaise('AIC is not meaningful for',self.kind)
        if self.edata is None:
            logger.TBRaise('AIC needs edata.')
        aic = AIC(self.xspline, self.yspline, self.cov, self)
        if self.deriv is not None:
            xd, dy, de, order = self.deriv
            for i in range(len(xd)):
                r    = (float(self(xd[i],der=order[i])) - dy[i])/de[i]
                aic += r**2 + np.log(2*np.pi*de[i]**2)
        return aic

    def __repr__(self) -> str:
        return "TBSpline"

    def __call__(self,x,der=0):
        return splev(x,self.tck,der=der)

    def get_knots(self):
        return self.tck[0]

    def get_coeffs(self):
        return self.tck[1]

    def get_order(self):
        return self.tck[2]

    def get_x(self):
        return self.xspline

    def get_y(self):
        return self.yspline

    def get_weights(self):
        """
        1/errors of the data, or None without edata.
        """
        if self.cov is None:
            return None
        return 1/np.sqrt(np.diag(self.cov))

    def n_deriv(self,x,n):
        """
        Take n derivatives of spline, evaluate at x

        Args:
            x (float)
            n (int)

        Returns:
            float: d^n S/ dx^n
        """
        return splev(x, self.tck, der=n)


#
# Kinds of spline. A kind holds the options of one kind of spline. Its fitter, made for given xdata and edata,
# does everything that does not depend on ydata once. fitter.fit(y) then returns a TBSpline. bootSpline calls
# fit once per bootstrap sample.
#


class SplineKind:

    """
    Base class of the kinds of spline. hasAIC is False for kinds without a meaningful AIC, i.e.
    interpolation and averages over models. natural is True for kinds with zero curvature at the endpoints.
    """

    hasAIC  = True
    natural = False

    def fitter(self, xdata, edata=None, deriv=None):
        logger.TBRaise('Not implemented for',self)


class Interpolating(SplineKind):

    """
    Natural cubic interpolating spline, with a knot at every datum. It passes through the data, so edata
    are not used in the fit.
    """

    hasAIC  = False
    natural = True

    def __repr__(self) -> str:
        return "Interpolating()"

    def fitter(self, xdata, edata=None, deriv=None):
        return _InterpolatingFitter(self, xdata, edata, deriv)


class Regression(SplineKind):

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

    def __init__(self, nknots=None, knots=None, rand=False, seed=None, natural=False, order=3):
        # Set the attributes first, so that __repr__ works in error messages.
        self.nknots  = nknots
        self.knots   = knots
        self.rand    = rand
        self.seed    = seed
        self.natural = natural
        self.order   = order
        if knots is None:
            self.knots = []
        self.knots = list(self.knots)
        if nknots is None:
            if knots is None:
                logger.TBRaise('Regression needs nknots or knots.')
            self.nknots = len(self.knots)
        checkType('int',nknots=self.nknots)
        checkType('int',order=order)
        if len(self.knots) > self.nknots:
            logger.TBRaise("len(knots)",len(self.knots),"exceeds nknots",self.nknots)
        if natural and (order != 3):
            logger.TBRaise("Natural splines have order=3 by definition.")
        if rand and (seed is None):
            self.seed = int(TBRNG().integers(2**32))

    def __repr__(self) -> str:
        return f"Regression(nknots={self.nknots}, natural={self.natural}, order={self.order})"

    def fitter(self, xdata, edata=None, deriv=None):
        return _RegressionFitter(self, xdata, edata, deriv)


class Smoothing(SplineKind):

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

    natural = True

    def __init__(self, lam='average', nlam=100):
        # Set the attributes first, so that __repr__ works in error messages.
        self.lam  = lam
        self.nlam = nlam
        if isinstance(lam,str):
            if lam != 'average':
                logger.TBRaise("lam must be a number or 'average'. lam =",lam)
            checkType('int',nlam=nlam)
            self.averaged = True
            self.hasAIC   = False
        else:
            checkType('real',lam=lam)
            if lam <= 0:
                logger.TBRaise('Need lam > 0. lam =',lam)
            self.averaged = False
            self.hasAIC   = True

    def __repr__(self) -> str:
        if self.averaged:
            return f"Smoothing(lam='average', nlam={self.nlam})"
        return f"Smoothing(lam={self.lam})"

    def fitter(self, xdata, edata=None, deriv=None):
        return _SmoothingFitter(self, xdata, edata, deriv)


class KnotAverage(SplineKind):

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

    hasAIC = False

    def __init__(self, kmax=None, ngrid=None, natural=False, nupdates=10000, nequi=1000, every=10, nchains=1, 
                 seed=DEFAULTSEED):
        # Set the attributes first, so that __repr__ works in error messages.
        self.kmax     = kmax
        self.ngrid    = ngrid
        self.natural  = natural
        self.nupdates = nupdates
        self.nequi    = nequi
        self.every    = every
        self.nchains  = nchains
        self.seed     = seed
        if kmax is not None:
            checkType('int',kmax=kmax)
            if kmax < 0:
                logger.TBRaise('Need kmax >= 0. kmax =',kmax)
        if ngrid is not None:
            checkType('int',ngrid=ngrid)
        checkType('int',nupdates=nupdates)
        checkType('int',nequi=nequi)
        checkType('int',every=every)
        checkType('int',nchains=nchains)
        checkType('int',seed=seed)
        if nequi < 0:
            logger.TBRaise('Need nequi >= 0. nequi =',nequi)
        if nupdates < 1:
            logger.TBRaise('Need nupdates >= 1. nupdates =',nupdates)

    def __repr__(self) -> str:
        return f"KnotAverage(kmax={self.kmax}, ngrid={self.ngrid}, natural={self.natural})"

    def fitter(self, xdata, edata=None, deriv=None):
        return _KnotAverageFitter(self, xdata, edata, deriv)


class _Fitter:

    """
    Base class of the fitters. Every fitter has fit(y, warn, dy), which returns a TBSpline; warn turns on
    warnings about the fit, where a fitter has any, and dy replaces the values of the derivative
    constraints, which bootSpline draws like data. Each subclass implements _fit(y, warn, dy).
    """

    # Highest derivative order a constraint may have
    maxDerivOrder = 3

    def __init__(self, kind, xdata, edata, deriv=None):
        self.kind  = kind
        self.x     = np.asarray(xdata,dtype=float)
        self.edata = edata
        self.xd, self.dyd, self.ed, self.od = _parseDeriv(self.x, deriv, edata, self.maxDerivOrder)

    def _derivOf(self, dy):
        """
        The derivative constraints with values dy, as stored in a TBSpline, or None without constraints.
        """
        if len(self.xd) == 0:
            return None
        return (self.xd, dy, self.ed, self.od)

    def _whitened(self, y, dy):
        """
        Whitened data M y, followed by the whitened values dy/de of the derivative constraints.
        """
        return np.r_[self.M @ y, np.asarray(dy,dtype=float)/self.ed]

    def fit(self, y, warn=False, dy=None):
        """
        Fit the data y, with the derivative constraints at the values dy, by default their given values.
        """
        if dy is None:
            dy = self.dyd
        return self._fit(y, warn, dy)

    def _fit(self, y, warn, dy):
        logger.TBRaise('Not implemented for',self.kind)


class _InterpolatingFitter(_Fitter):

    def __init__(self, kind, xdata, edata, deriv=None):
        super().__init__(kind, xdata, edata, deriv)
        if len(self.xd) > 0:
            logger.TBRaise('Interpolating passes through the data, so it cannot take derivative constraints.')

    def _fit(self, y, warn, dy):
        spl = make_interp_spline(self.x, y, k=3, bc_type='natural')
        return TBSpline(self.x, y, (spl.t, spl.c, 3), self.kind, edata=self.edata)


def _knotVector(xmin, xmax, interior, order=3):
    """
    Full knot vector as in splrep: each endpoint repeated order+1 times around the interior knots.
    """
    return np.r_[[xmin]*(order+1), interior, [xmax]*(order+1)]


def _lsqSystem(x, M, t, order, natural, xd, od, ed):
    """
    Set up the least-squares problem |X d - b|^2 of a regression spline with full knot vector t, where
    X = M B N, b = M y, and c = N d. For a natural spline, N imposes S''(a) = S''(b) = 0, see _naturalBasis;
    otherwise N = 1. Each derivative constraint adds the row
    B^(od)(xd) N / ed to X. The solution d solves the normal equations A d = X^T b with A = X^T X. Returns
    N, X, A, and whether the fit cannot be solved: A is singular when the Schoenberg-Whitney conditions fail,
    and too ill-conditioned when derivative constraints are much more precise than the data.
    """
    if natural:
        N = _naturalBasis(t)
    else:
        N = np.eye(len(t) - order - 1)
    X = _fitRows(x, M, t, order, N, xd, od, ed)
    A = X.T @ X
    # With fewer rows than columns, A is singular, and its condition number is infinite.
    return N, X, A, np.linalg.cond(A) > _CONDMAX


class _RegressionFitter(_Fitter):

    def __init__(self, kind, xdata, edata, deriv=None):
        """
        Place the knots and set up the least-squares problem, see _lsqSystem.
        """
        self.maxDerivOrder = kind.order
        super().__init__(kind, xdata, edata, deriv)
        order      = kind.order
        xmin, xmax = np.min(self.x), np.max(self.x)
        nauto = kind.nknots - len(kind.knots)
        if nauto > 0:
            if kind.rand:
                knots = _random_knots(self.x,nauto,SEED=kind.seed)
            else:
                knots = _even_knots(self.x,nauto)
        else:
            knots = []
        knots = sorted(knots + kind.knots)
        if len(knots)>0 and knots[0]<xmin:
            logger.TBRaise("You can't put a knot to the left of the x-data. knots, min(xdata) = ",knots,xmin)
        if len(knots)>0 and knots[-1]>xmax:
            logger.TBRaise("You can't put a knot to the right of the x-data. knots, max(xdata) = ",knots,xmax)
        self.t = _knotVector(xmin, xmax, knots, order)
        self.M = whiteningMatrix(parseEdata(edata, len(self.x)))
        self.N, self.X, self.A, singular = _lsqSystem(self.x, self.M, self.t, order, kind.natural, self.xd, self.od,
                                                      self.ed)
        if singular:
            logger.TBRaise('The fit cannot be solved. Each interval between knots needs enough data (Schoenberg-Whitney),',
                           'and derivative constraints must not be much more precise than the data.')

    def _fit(self, y, warn, dy):
        d = np.linalg.solve(self.A, self.X.T @ self._whitened(y, dy))
        tck = (self.t, self.N @ d, self.kind.order)
        return TBSpline(self.x, y, tck, self.kind, edata=self.edata, deriv=self._derivOf(dy))


class _SmoothingFitter(_Fitter):

    def __init__(self, kind, xdata, edata, deriv=None):
        """
        Set up the system of _smoothingSystem, and _smoothingAtLam for each lam: the given one, or the grid for
        lam='average'. Everything except the weights is linear in y.
        """
        super().__init__(kind, xdata, edata, deriv)
        self.t, self.M, self.X, self.Omega = _smoothingSystem(self.x, edata, self.xd, self.od, self.ed)
        if kind.averaged:
            if edata is None:
                logger.TBRaise("Averaging over lam needs edata, since the evidence depends on the errors.")
            self.lams = _lamGrid(self.X, self.Omega, len(np.unique(self.x)), len(self.xd), kind.nlam)
        else:
            self.lams = np.array([kind.lam])
        A, Ainv, peffs, occam = [], [], [], []
        for lam in self.lams:
            Ak, Ainvk, peffk, occamk = _smoothingAtLam(self.X, self.Omega, lam)
            A.append(Ak)
            Ainv.append(Ainvk)
            peffs.append(peffk)
            occam.append(occamk)
        self.A     = np.array(A)
        self.Ainv  = np.array(Ainv)
        self.peffs = np.array(peffs)
        self.occam = np.array(occam)

    def _coeffs(self, y, warn, dy):
        """
        Weights, log evidences, and coefficients at each lam.
        """
        w, logEv, c = _averageWeights(self._whitened(y, dy), self.X, self.Omega, self.lams, self.A, self.occam)
        if warn and self.kind.averaged and ((w[0] > 0.01) or (w[-1] > 0.01)):
            logger.warn('Weights at the ends of the lam grid are',w[0],'and',w[-1],'which may make the result',
                        'depend on the grid. Large weight at large lam means the data are consistent with a straight line.')
        return w, logEv, c

    def _fit(self, y, warn, dy):
        """
        All splines share the knots, so their weighted average is again a spline, whose coefficients are the
        weighted average of the coefficients.
        """
        w, logEv, c = self._coeffs(y, warn, dy)
        tck = (self.t, w @ c, 3)
        if self.kind.averaged:
            return TBSpline(self.x, y, tck, self.kind, edata=self.edata, peff=w @ self.peffs,
                            deriv=self._derivOf(dy), lam='average')
        return TBSpline(self.x, y, tck, self.kind, edata=self.edata, peff=self.peffs[0], deriv=self._derivOf(dy),
                        lam=self.lams[0], logEvidence=logEv[0])


def _naturalBasis(t):
    """
    Matrix N with c = N d for the cubic splines with S''(a) = S''(b) = 0 and full knot vector t. The free
    coefficients d are c_1,...,c_{nB-2}; the conditions fix c_0 and c_{nB-1}.
    """
    nB = len(t) - 4
    N  = np.zeros((nB, nB-2))
    for j in range(nB-2):
        N[j+1,j] = 1.
    # S''(a) = 0 means equal slopes of the control polygon on its first two legs, likewise at b.
    a = (t[4]-t[1])/(t[5]-t[2])
    N[0,:] = N[1,:] - a*(N[2,:]-N[1,:])
    b = (t[nB+2]-t[nB-1])/(t[nB+1]-t[nB-2])
    N[nB-1,:] = N[nB-2,:] + b*(N[nB-2,:]-N[nB-3,:])
    return N


# The Metropolis-Hastings chain of KnotAverage. A configuration is the sorted list of grid indices of its knots.


def _knotValid(idx, lt, le, nu):
    """
    Is there a unique datum strictly inside each interval of the configuration idx?
    """
    k = len(idx)
    if k == 0:
        return True
    if lt[idx[0]] < 2:
        return False
    if nu - le[idx[k-1]] < 2:
        return False
    for i in range(k-1):
        if lt[idx[i+1]] - le[idx[i]] < 1:
            return False
    return True


def _knotFit(x, M, My, grid, idx, natural, xd, od, ed):
    """
    Least-squares fit of the cubic regression spline with the configuration idx to the whitened data My, see
    _lsqSystem. The last len(xd) entries of My are the whitened derivative constraints. Returns whether the fit
    is determined, its AIC, and its B-spline coefficients.
    """
    t = _knotVector(x.min(), x.max(), grid[idx])
    N, X, A, singular = _lsqSystem(x, M, t, 3, natural, xd, od, ed)
    if singular:
        return False, 0., np.zeros(len(idx)+4)
    d    = np.linalg.solve(A, X.T @ My)
    r    = My - X @ d
    chi2 = r @ r
    return True, chi2 + 2*X.shape[1], N @ d


def _freePoint(idx, m):
    """
    The m-th grid index, counting from 0, that is not in the configuration idx.
    """
    g = m
    for i in idx:
        if i <= g:
            g += 1
        else:
            break
    return g


def _knotChain(x, M, My, grid, lt, le, nu, kmax, natural, idx0, u, nequi, every, wmove, pbirth, pdeath, xd, od, ed):
    """
    Metropolis-Hastings chain over knot configurations, starting from the configuration idx0, with u[it,:] the
    6 uniform random numbers of update it. The first nequi updates are for thermalization; after them, every
    every-th configuration is saved. xd, od, ed are the derivative constraints, see _knotFit. Returns a list
    with the configuration, AIC, and coefficients of every saved configuration.
    """
    G   = len(grid)
    idx = list(idx0)
    ok, aic, c = _knotFit(x, M, My, grid, idx, natural, xd, od, ed)
    saved = []
    for it in range(len(u)):
        k = len(idx)
        if k < kmax:
            b = pbirth
        else:
            b = 0.
        if k > 0:
            d = pdeath
        else:
            d = 0.
        new = None
        if u[it,0] < b:
            # Birth at a random free grid point. Proposal ratio d_{k+1}/b_k; the grid factors cancel
            # against p(t|k), see the notes.
            g        = _freePoint(idx, min(int(u[it,1]*(G-k)), G-k-1))
            new      = sorted(idx + [g])
            logRatio = np.log(pdeath/b)
        elif u[it,0] < b + d:
            j        = min(int(u[it,1]*k), k-1)
            new      = idx[:j] + idx[j+1:]
            logRatio = np.log(pbirth/d)
        elif k > 0:
            j        = min(int(u[it,1]*k), k-1)
            logRatio = 0.
            if u[it,2] < 0.5:
                # Relocate knot j to a random grid point not occupied by any knot.
                g = _freePoint(idx, min(int(u[it,3]*(G-k)), G-k-1))
            else:
                # Shift knot j by up to wmove grid points.
                step = 1 + min(int(u[it,3]*wmove), wmove-1)
                if u[it,4] < 0.5:
                    step = -step
                g = idx[j] + step
                if (g < 0) or (g >= G) or (g in idx):
                    g = None
            if g is not None:
                new = sorted(idx[:j] + idx[j+1:] + [g])
        if (new is not None) and _knotValid(new, lt, le, nu):
            oknew, aicnew, cnew = _knotFit(x, M, My, grid, new, natural, xd, od, ed)
            if oknew:
                logA = -0.5*(aicnew - aic) + logRatio
                if (logA >= 0) or ((logA > -700) and (u[it,5] < np.exp(logA))):
                    idx = new
                    aic = aicnew
                    c   = cnew
        if (it >= nequi) and ((it - nequi) % every == 0):
            saved.append((tuple(idx), aic, c))
    return saved


class _KnotAverageFitter(_Fitter):

    # Probability to propose a birth (death) when k < kmax (k > 0); the rest are moves.
    PBIRTH = 1/3
    PDEATH = 1/3
    # Warn if the Geweke z of the AIC trace exceeds DRIFTZ, or if a chain saved fewer than MINSAVED configurations.
    DRIFTZ   = 3
    MINSAVED = 50

    def __init__(self, kind, xdata, edata, deriv=None):
        """
        Set up the grid of candidate knots, the largest number of knots, and the basis with a knot at every 
        grid point, in which the averaged spline is expressed.
        """
        if edata is None:
            logger.TBRaise('KnotAverage needs edata, since the AIC depends on the errors.')
        super().__init__(kind, xdata, edata, deriv)
        self.M     = whiteningMatrix(parseEdata(edata, len(self.x)))
        u          = np.unique(self.x)
        nu         = len(u)
        self.xmin, self.xmax = u[0], u[-1]
        if kind.ngrid is None:
            self.ngrid = max(100, int(np.ceil(4*(self.xmax-self.xmin)/np.min(np.diff(u)))))
        else:
            self.ngrid = kind.ngrid
        self.grid = np.linspace(self.xmin, self.xmax, self.ngrid+2)[1:-1]
        # Number of unique x strictly below, and below or at, each grid point
        self.lt = np.searchsorted(u, self.grid, side='left')
        self.le = np.searchsorted(u, self.grid, side='right')
        self.nu = nu
        # A configuration needs a datum strictly inside each of its k+1 intervals, k <= nu-3, and at most 
        # as many parameters as unique x.
        if kind.natural:
            ceiling = nu-3
        else:
            ceiling = nu-4
        if ceiling < 0:
            logger.TBRaise('Too few unique x for a cubic spline. Unique x:',nu)
        # Every knot needs its own grid point.
        ceiling = min(ceiling, self.ngrid)
        if kind.kmax is None:
            self.kmax = ceiling
        else:
            if kind.kmax > ceiling:
                logger.TBRaise('kmax =',kind.kmax,'exceeds the largest number of knots the data allow,',ceiling)
            self.kmax = kind.kmax
        # Basis with a knot at every grid point. Every configuration is a spline in this basis, so its 
        # coefficients follow from its values at the Greville abscissae z of this basis.
        self.tFull = _knotVector(self.xmin, self.xmax, self.grid)
        # Clip, since rounding can put the end abscissae just outside [xmin, xmax], e.g. (3*1.6)/3 > 1.6.
        self.z     = np.clip((self.tFull[1:-3] + self.tFull[2:-2] + self.tFull[3:-1])/3, self.xmin, self.xmax)
        self.Pz    = np.linalg.inv(BSpline.design_matrix(self.z, self.tFull, 3).toarray())
        self._cacheKey = None
        self._cacheRun = None

    def _system(self, idx, My):
        """
        Fit of the configuration with grid indices idx to the whitened data My, see _whitened, by the same 
        _knotValid and _knotFit as the chain. Returns its AIC, knot vector, coefficients, and number of parameters, 
        or None if the configuration is not allowed.
        """
        idx = sorted(int(i) for i in idx)
        if not _knotValid(idx, self.lt, self.le, self.nu):
            return None
        ok, aic, c = _knotFit(self.x, self.M, My, self.grid, idx, self.kind.natural, self.xd, self.od, self.ed)
        if not ok:
            return None
        t = _knotVector(self.xmin, self.xmax, self.grid[idx])
        return {'AIC':aic, 't':t, 'c':c, 'p':_nparams(len(idx), self.kind.natural)}

    def _chain(self, My, rng, saved, trace):
        """
        One Metropolis-Hastings chain for the whitened data My, run by _knotChain. Adds the saved configurations 
        to saved, a dict from the tuple of grid indices to [count, fit], and the AIC of every saved
        configuration, in order, to trace.
        """
        kind = self.kind
        # Start from up to 2 knots, evenly spaced in data index, moved to the nearest grid points.
        idx0 = []
        if min(2, self.kmax) > 0:
            start = np.unique(np.argmin(np.abs(np.subtract.outer(_even_knots(self.x,min(2,self.kmax)), self.grid)),axis=1))
            if self._system(start, My) is not None:
                idx0 = start.tolist()
        # 6 uniform random numbers per update: move type, which knot or grid point, relocate or shift, where to,
        # direction of the shift, and acceptance.
        u = rng.random((kind.nequi + kind.nupdates, 6))
        # Local moves shift a knot by up to wmove grid points.
        wmove = max(1, self.ngrid//20)
        chain = _knotChain(self.x, self.M, My, self.grid, self.lt, self.le, self.nu, self.kmax, kind.natural, idx0, u,
                           kind.nequi, kind.every, wmove, self.PBIRTH, self.PDEATH, self.xd, self.od, self.ed)
        for key, aic, c in chain:
            trace.append(aic)
            if key in saved:
                saved[key][0] += 1
            else:
                t = _knotVector(self.xmin, self.xmax, self.grid[list(key)])
                saved[key] = [1, {'AIC':aic, 't':t, 'c':c, 'p':_nparams(len(key), kind.natural)}]

    def _run(self, y, dy):
        """
        Run the chains for the data y and derivative constraints dy and pool their saved configurations. The 
        result for the last y and dy is cached, since details needs it again.
        """
        yall = np.r_[np.asarray(y,dtype=float), np.asarray(dy,dtype=float)]
        key  = yall.tobytes()
        if key == self._cacheKey:
            return self._cacheRun
        My    = self._whitened(y, dy)
        saved  = {}
        traces = []
        # Seed each chain from the seed of the kind, the chain index, and the bits of the data.
        words = np.frombuffer(np.ascontiguousarray(yall,dtype=np.float64).tobytes(),dtype=np.uint32)
        words = words.astype(np.int64).tolist()
        for ichain in range(self.kind.nchains):
            rng   = TBRNG([self.kind.seed, ichain] + words)
            trace = []
            self._chain(My, rng, saved, trace)
            traces.append(np.array(trace))
        counts, fits = [], []
        for cfg in saved:
            counts.append(saved[cfg][0])
            fits.append(saved[cfg][1])
        counts = np.array(counts)
        # How often each number of knots 0,...,kmax was saved. The knot vector of a fit with k interior knots has
        # k+8 entries.
        kHist  = np.zeros(self.kmax+1,dtype=int)
        for i in range(len(fits)):
            kHist[len(fits[i]['t'])-8] += counts[i]
        run = {'counts':counts, 'fits':fits, 'traces':traces, 'kHistogram':kHist}
        self._cacheKey = key
        self._cacheRun = run
        return run

    def _fit(self, y, warn, dy):
        """
        The chain visits configurations with frequency proportional to their weight, so the average is the
        plain average over the saved configurations, weighted by how often each was saved.
        """
        run    = self._run(y, dy)
        counts = run['counts']
        sz     = np.zeros(len(self.z))
        peff   = 0.
        for i in range(len(counts)):
            fit   = run['fits'][i]
            sz   += counts[i]*(BSpline.design_matrix(self.z, fit['t'], 3).toarray() @ fit['c'])
            peff += counts[i]*fit['p']
        ntot  = np.sum(counts)
        sz   /= ntot
        peff /= ntot
        atMax = run['kHistogram'][self.kmax]/ntot
        if warn and (atMax > 0.01):
            logger.warn('The chain spends',atMax,'of its time at kmax =',self.kmax,'knots. The data',
                        'may want more knots than allowed, or the configurations are fitting noise.')
        if warn:
            self._checkDrift(run['traces'])
        tck = (self.tFull, self.Pz @ sz, 3)
        return TBSpline(self.x, y, tck, self.kind, edata=self.edata, peff=peff, deriv=self._derivOf(dy))

    def _checkDrift(self, traces):
        """
        Geweke test of each chain, see gewekeZ: compare the mean AIC of the first 10% of its saved configurations
        with the mean AIC of the last 50%. A chain still drifting after nequi gives a large |z|.
        """
        for ichain in range(len(traces)):
            trace = traces[ichain]
            if len(trace) < self.MINSAVED:
                logger.warn('Chain',ichain,'saved only',len(trace),'configurations, too few to check whether it',
                            'thermalized. Increase nupdates.')
                continue
            z = gewekeZ(np.asarray(trace,dtype=float))
            if np.abs(z) > self.DRIFTZ:
                logger.warn('Chain',ichain,'may not be thermalized: the mean AIC of the first 10% of its saved',
                            'configurations differs from that of the last 50% by z =',round(z,1),'standard errors.',
                            'Increase nequi and nupdates.')

    def details(self, y, xspl, dy):
        """
        The spline and weight of each saved configuration for the data y and derivative constraints dy, and the
        weighted spreads over configurations of the spline and its observables on xspl, see _modelSpreads.
        """
        run    = self._run(y, dy)
        counts = run['counts']
        w      = counts/np.sum(counts)
        splines = []
        for i in range(len(counts)):
            fit     = run['fits'][i]
            cfgKind = Regression(knots=list(fit['t'][4:-4]), natural=self.kind.natural)
            splines.append(TBSpline(self.x, y, (fit['t'], fit['c'], 3), cfgKind, edata=self.edata,
                                    deriv=self._derivOf(dy)))
        det = {}
        det['splines'] = splines
        det['weights'] = w
        det['spread']  = _modelSpreads(splines, w, xspl)
        return det


def getSpline(xdata, ydata, edata=None, kind=Interpolating(), deriv=None):
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
    checkEqualLengths(xdata,ydata)
    checkEdata(edata, len(xdata))
    return kind.fitter(xdata, edata, deriv).fit(np.asarray(ydata,dtype=float), warn=True)


def getSplineErr(xdata, xspline, ydata, ydatae, kind=Interpolating()):
    """
    Fit unweighted splines of the given kind to ydata-ydatae and ydata+ydatae, evaluate them at xspline,
    and return their midpoint and half-difference. This gives a smooth band, but it is not a statistical
    error propagation; for that use bootSpline.
    """
    if np.ndim(ydatae) != 1:
        logger.TBRaise('getSplineErr needs a vector of errors ydatae.')
    spline_lower  = getSpline(xdata, ydata - ydatae, kind=kind)(xspline)
    spline_upper  = getSpline(xdata, ydata + ydatae, kind=kind)(xspline)
    spline_center = (spline_lower+spline_upper)/2
    spline_err    = (spline_upper-spline_lower)/2
    return spline_center, spline_err


def _steepestInflection(spline, xspl):
    """
    The inflection point (x,y) of spline on xspl where |spline'| is largest. Inflection points are where the
    second derivative changes sign; each root is found between the two grid points that bracket it. Taking
    the steepest one ignores the small wiggles a smoothed spline can have where the data are flat. A zero of
    the second derivative at an end of xspl, as for a natural spline, does not count, so the first and last
    grid intervals are skipped; numerically that zero is not exact. Returns (nan,nan) if there is no
    inflection point.
    """
    d2     = spline(xspl,der=2)
    change = np.where(d2[1:-2]*d2[2:-1] < 0)[0] + 1
    if len(change) == 0:
        return np.nan, np.nan
    xbest, slopebest = np.nan, -np.inf
    for i in change:
        x     = brentq(lambda t: spline(t,der=2), xspl[i], xspl[i+1])
        slope = abs(float(spline(x,der=1)))
        if slope > slopebest:
            xbest, slopebest = x, slope
    return xbest, float(spline(xbest))


def _observables(spline, xspl):
    """
    The observables of a spline that bootSpline and posteriorSpline report, keyed by their names in the result:
    the spline on xspl ('yspl'), the x on xspl where it is largest and its value there ('xmax', 'ymax'), and
    its steepest inflection point ('xinflect', 'yinflect'), which is nan if there is none.
    """
    ys   = spline(xspl)
    imax = np.argmax(ys)
    xinf, yinf = _steepestInflection(spline, xspl)
    return {'yspl':ys, 'xmax':xspl[imax], 'ymax':ys[imax], 'xinflect':xinf, 'yinflect':yinf}


def _values(obs, name):
    """
    The values of the observable name in the list obs of observables, see _observables.
    """
    vals = []
    for o in obs:
        vals.append(o[name])
    return vals


def _defined(vals, name, what):
    """
    The entries of vals without nan. Warns how many of the what (e.g. samples) were dropped because the
    observable name is undefined for them.
    """
    keep = []
    for val in vals:
        if not np.any(np.isnan(val)):
            keep.append(val)
    if len(keep) < len(vals):
        logger.warn(f'Dropped {len(vals)-len(keep)} of {len(vals)} {what} where {name} is undefined.')
    return np.array(keep)


def _spreadOf(vals, w):
    """
    Weighted spread over models, with weights w, of their values vals of a quantity, which may be a number or
    an array. Models where the quantity has a nan are dropped, and the remaining weights are renormalized.
    Returns nan if it has a nan for every model.
    """
    keep, wkeep = [], []
    for i in range(len(vals)):
        val = np.asarray(vals[i],dtype=float)
        if not np.any(np.isnan(val)):
            keep.append(val)
            wkeep.append(w[i])
    if len(wkeep) == 0:
        return np.nan
    wkeep = np.array(wkeep)
    if np.sum(wkeep) < 0.99:
        logger.warn(f'Models where the quantity is undefined carry weight {1-np.sum(wkeep):.3f}; they are dropped from its systematic error.')
    return weightedSpread(wkeep, np.array(keep))


def _modelSpreads(splines, w, xspl):
    """
    Systematic errors of an average over models, e.g. knot configurations or values of lam: the weighted
    spreads over the models, with weights w, of each observable of their splines, see _observables. Models
    where an observable is undefined are dropped from its spread, see _spreadOf.
    """
    obs = []
    for spline in splines:
        obs.append(_observables(spline, xspl))
    sp = {}
    for name in obs[0]:
        sp[name] = _spreadOf(_values(obs, name), w)
    return sp


def modelSpread(res, f):
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
    if 'splinesModel' not in res:
        return 0.
    vals = []
    for spline in res['splinesModel']:
        vals.append(f(spline))
    return _spreadOf(vals, res['weightsModel'])


def _kindSample(sample, fitter):
    """
    Fit one bootstrap sample for bootSpline: ydata, followed by the values of the derivative constraints.
    The bootstrap needs numbers, so pack what is needed to rebuild the spline into one flat array: the 
    sample, the coefficients, and the number of parameters. bootSpline unpacks it. This lives at module 
    level so that it can be pickled for nproc>1.
    """
    yBS = sample[0]
    n   = len(fitter.x)
    spl = fitter.fit(yBS[:n], dy=yBS[n:])
    return np.r_[yBS, spl.get_coeffs(), spl.get_nparams()]


def _setErrors(res, name, ename, f0, fBS, sys):
    """
    Put the bias-corrected central value of an observable into res[name], and its statistical error over
    the bootstrap samples fBS, its systematic error sys, and their sum in quadrature into res[ename+'_stat'],
    res[ename+'_sys'], and res[ename].
    """
    stat = dev_by_dist(fBS)
    res[name]          = biasCorrect(f0, fBS, average='median')
    res[ename+'_stat'] = stat
    res[ename+'_sys']  = sys
    res[ename]         = quadrature(np.array([stat, sys]),axis=0)


def bootSpline(xdata, ydata, edata, kind=Interpolating(), numb_samples=300, nsupport=301, seed=DEFAULTSEED, nproc=1, 
               deriv=None) -> dict:
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
        For KnotAverage also:
            splinesModel, weightsModel: the spline of each model (knot configuration) for the original
                data and its weight. Use modelSpread for the systematic error of your own observables.
    """
    checkType("int",seed=seed)
    checkEqualLengths(xdata,ydata)
    if edata is None:
        logger.TBRaise('bootSpline draws its samples from edata, so it needs edata. For a spline without errors, use getSpline.')
    checkEdata(edata, len(xdata))
    y = np.asarray(ydata,dtype=float)
    cov = parseEdata(edata)
    std = np.sqrt(np.diag(cov))
    if isMatrix(np.asarray(edata)):
        useCovariance = True
        Covariance    = [cov]
    else:
        useCovariance = False
        Covariance    = None
    xspl  = np.linspace(np.min(xdata),np.max(xdata),nsupport)
    ndata = len(y)

    fitter  = kind.fitter(xdata, edata, deriv)
    spline0 = fitter.fit(y, warn=True)
    t, _, k = spline0.tck

    # The derivative constraints are drawn like data, so append them to the data.
    nsample = ndata + len(fitter.xd)
    ysample = np.r_[y, fitter.dyd]
    std     = np.r_[std, fitter.ed]
    if useCovariance:
        Covariance = [block_diag(cov, np.diag(fitter.ed**2))]

    samples, _, _ = bootstr_from_gauss(_kindSample,data=[ysample],data_std_dev=[std],
                                       numb_samples=numb_samples,return_sample=True,seed=seed,
                                       args={'fitter':fitter},nproc=nproc,
                                       useCovariance=useCovariance,Covariance=Covariance)
    samples = np.asarray(samples)

    # Rebuild the bootstrap splines without refitting.
    splines = []
    for sample in samples:
        # Unpack what _kindSample packed: the sample (ydata and derivative values), coefficients, number of parameters
        yBS     = sample[:nsample]
        c       = sample[nsample:-1]
        nparams = sample[-1]
        splines.append(TBSpline(xdata,yBS[:ndata],(t,c,k),kind,edata=edata,peff=nparams,
                                deriv=fitter._derivOf(yBS[ndata:]),lam=spline0.lam))

    # Systematic errors from the spread over models for the original data. Kinds that are not averages have none.
    if isinstance(kind,KnotAverage):
        det    = fitter.details(y, xspl, fitter.dyd)
        spread = det['spread']
    else:
        det    = None
        spread = None

    res = {}
    res['splineBS']   = splines # in case you want spline functions at bootstrap level
    res['xspl']       = xspl
    res['splineMean'] = spline0
    # Central values are 2 f(original data) - median over samples. The statistical errors are dev_by_dist over 
    # samples, and the totals add the systematic ones in quadrature. spline0 is the fit to the original data,
    # which is what the bootstrap's zero-noise sample would give. Because of the bias correction, xmax can lie
    # between grid points. Samples where an observable is undefined, e.g. without an inflection point, are
    # dropped from its error.
    obs0  = _observables(spline0, xspl)
    obsBS = []
    for spline in splines:
        obsBS.append(_observables(spline, xspl))
    for name in obs0:
        fBS = _defined(_values(obsBS, name), name, 'bootstrap samples')
        if np.any(np.isnan(obs0[name])) or len(fBS) == 0:
            for key in [name, name+'e', name+'e_stat', name+'e_sys']:
                res[key] = np.nan
            continue
        if spread is None:
            sys = np.zeros(np.shape(obs0[name]))
        else:
            sys = spread[name]
        _setErrors(res, name, name+'e', obs0[name], fBS, sys)
    if kind.hasAIC:
        res['AIC'] = spline0.get_AIC()
    if isinstance(kind,KnotAverage):
        res['splinesModel'] = det['splines']
        res['weightsModel'] = det['weights']
    return res


def posteriorSpline(xdata, ydata, edata, kind=Smoothing(), numb_samples=300, nsupport=301, seed=DEFAULTSEED,
                    deriv=None) -> dict:
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
    checkType("int",seed=seed)
    checkType("int",numb_samples=numb_samples)
    checkEqualLengths(xdata,ydata)
    checkEdata(edata, len(xdata))
    if not isinstance(kind,Smoothing):
        logger.TBRaise('posteriorSpline needs a Smoothing kind. Got',kind)
    y       = np.asarray(ydata,dtype=float)
    xspl    = np.linspace(np.min(xdata),np.max(xdata),nsupport)
    fitter  = kind.fitter(xdata, edata, deriv)
    spline0 = fitter.fit(y, warn=True)
    w, _, c = fitter._coeffs(y, False, fitter.dyd)
    nlam    = len(fitter.lams)
    nB      = fitter.Ainv.shape[2]
    # The posterior mean at each lam, and the systematic errors from its spread over lam
    splinesLam = []
    for k in range(nlam):
        splinesLam.append(TBSpline(xdata,y,(fitter.t,c[k],3),kind,edata=edata,peff=fitter.peffs[k],
                                   deriv=spline0.deriv,lam=fitter.lams[k]))
    spread = _modelSpreads(splinesLam, w, xspl)

    # At fixed lam, the posterior covariance of the coefficients is A^-1, so the posterior variance of the
    # spline on xspl is the diagonal of Bspl A^-1 Bspl^T. Average it over lam.
    Bspl  = BSpline.design_matrix(xspl, fitter.t, 3).toarray()
    ystat = np.zeros(len(xspl))
    for k in range(nlam):
        ystat += w[k]*np.sum((Bspl @ fitter.Ainv[k])*Bspl,axis=1)
    ystat = np.sqrt(ystat)

    # Draw lam_k with probability w_k, then the coefficients from the posterior at lam_k
    rng = TBRNG(seed)
    ks  = rng.choice(nlam, size=numb_samples, p=w)
    cs  = np.zeros((numb_samples, nB))
    for k in range(nlam):
        drawn = np.where(ks == k)[0]
        if len(drawn) > 0:
            cs[drawn] = rng.multivariate_normal(c[k], fitter.Ainv[k], size=len(drawn))
    splines, obs = [], []
    for i in range(numb_samples):
        spline = TBSpline(xdata,y,(fitter.t,cs[i],3),kind,edata=edata,peff=spline0.peff,deriv=spline0.deriv,
                          lam=spline0.lam)
        splines.append(spline)
        obs.append(_observables(spline, xspl))

    res = {}
    res['xspl']          = xspl
    res['yspl']          = spline0(xspl)
    res['ysple_stat']    = ystat
    res['ysple_sys']     = spread['yspl']
    res['ysple']         = quadrature(np.array([ystat, spread['yspl']]),axis=0)
    res['splineMean']    = spline0
    res['splineSamples'] = splines
    for name in ['xmax','ymax','xinflect','yinflect']:
        vals = _defined(_values(obs, name), name, 'posterior samples')
        if len(vals) == 0:
            for key in [name, name+'e', name+'e_stat', name+'e_sys']:
                res[key] = np.nan
            continue
        res[name]          = std_median(vals)
        res[name+'e']      = dev_by_dist(vals)
        res[name+'e_sys']  = spread[name]
        res[name+'e_stat'] = np.sqrt(max(res[name+'e']**2 - spread[name]**2, 0))
    res['lams']          = fitter.lams
    res['peffs']         = fitter.peffs
    res['weights']       = w
    return res
