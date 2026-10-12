#
# spline.py
#
# D. Clarke, Claude Code
#
# Splines with control over knot placement, endpoints, and smoothing, and their error bands. What kind of
# spline you get is set by a kind object: Interpolating, Regression, Smoothing, or KnotAverage. Pass it to
# getSpline for a single fit, to bootSpline for a bootstrap error band, or to posteriorSpline for a Bayesian
# one.
#

import numpy as np
from scipy.interpolate import BSpline, make_interp_spline, splev
from scipy.linalg import null_space, solve_triangular
from scipy.optimize import brentq
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import AIC, whiteningMatrix, getModelWeights
from latqcdtools.base.check import checkType, checkEqualLengths
from latqcdtools.base.initialize import TBRNG, DEFAULTSEED
from latqcdtools.base.speedify import compile
from latqcdtools.math.math import isMatrix, isVector, quadrature
from latqcdtools.statistics.statistics import std_median, dev_by_dist
from latqcdtools.statistics.bootstr import bootstr_from_gauss


def _even_knots(xdata, nknots):
    """
    Return a list of nknots knots, evenly spaced in data index (not in x), so that each interval holds
    roughly the same number of data. Each knot is the midpoint of two neighboring unique data, or a datum
    itself when the index lands exactly on one.
    """
    if len(xdata)<nknots:
        logger.TBRaise('number of data < number of knots')
    flat_xdata = np.sort(np.asarray(xdata))
    # remove duplicate x values
    flat_xdata = np.unique(flat_xdata)
    jump_step = (len(flat_xdata) - 1) / (nknots + 1)
    knots = []
    for i in range(1, nknots + 1):
        x_lower = flat_xdata[int(np.floor(i * jump_step))]
        x_upper = flat_xdata[int(np.ceil(i * jump_step))]
        knots.append((x_lower + x_upper) / 2)
    return knots


def _random_knots(xdata, nknots, randomization_factor=1, SEED=None):
    """
    Return a list of nknots randomly placed knots. Draws a random subset of xdata and applies _even_knots
    to it. randomization_factor=1 draws the smallest subset (most random); 0 uses all the data.
    """
    rng = TBRNG(SEED)
    flat_xdata = np.sort(np.asarray(xdata))
    nsample = int(nknots+1+(1-randomization_factor)*(len(flat_xdata)-nknots))
    sample_xdata = rng.choice(flat_xdata,nsample,replace=False)
    # Retry if too many data points are removed by np.unique
    while len(np.unique(sample_xdata)) < nknots + 1:
        sample_xdata = rng.choice(flat_xdata,nsample,replace=False)
    return _even_knots(sample_xdata, nknots)


def _parseEdata(edata):
    """
    Interpret edata like the Fitter does: a vector holds errors and a matrix is the covariance matrix.
    Returns the weights 1/sigma and the covariance matrix, which is None unless edata is a matrix.
    """
    if edata is None:
        return None, None
    edata = np.asarray(edata,dtype=float)
    if isVector(edata):
        return 1/edata, None
    elif isMatrix(edata):
        return 1/np.sqrt(np.diag(edata)), edata
    else:
        logger.TBRaise('edata must be a vector of errors or a covariance matrix. Got ndim =',edata.ndim)


def _whitening(edata, ndata):
    """
    Whitening matrix M with M^T M = W, so that a chi^2 (y - S)^T W (y - S) becomes |M y - M S|^2. For
    errors, M = diag(1/edata). For a covariance matrix, M is the whiteningMatrix of the toolbox. Without
    edata, M = 1.
    """
    weights, cov = _parseEdata(edata)
    if cov is not None:
        return whiteningMatrix(cov)
    elif weights is None:
        return np.eye(ndata)
    else:
        return np.diag(weights)


def _smoothingSystem(xdata, edata, lam):
    """
    Set up and solve by QR the least-squares problem of a cubic smoothing spline, which minimizes
    (y - S)^T W (y - S) + lam int S''^2. The minimizer over all twice-differentiable functions is a natural
    cubic spline with a knot at every unique x, so we work in that B-spline basis. With S = B c and
    int S''^2 = |L c|^2, the problem is |X c - b|^2 with X = (M B, sqrt(lam) L) and b = (M y, 0) stacked,
    where M is the whitening matrix. QR, X = Q R, stays accurate for large lam, unlike the normal equations.
    None of this depends on y.

    Returns:
        t: full knot vector
        M: whitening matrix
        Qd: rows of Q belonging to the data, so that c = R^-1 Qd^T M y
        R: triangular factor
        occam: the part of the log evidence that does not depend on y, see _logEvidence
    """
    x = np.asarray(xdata,dtype=float)
    u = np.unique(x)
    t = np.r_[[u[0]]*4, u[1:-1], [u[-1]]*4]
    nB = len(t) - 4
    B  = BSpline.design_matrix(x, t, 3).toarray()
    M  = _whitening(edata, len(x))
    # B_j'' is linear between neighboring knots, so B_i'' B_j'' is quadratic there, and 2-point
    # Gauss-Legendre quadrature on each interval gives int S''^2 = sum_q wq S''(xq)^2 exactly.
    half = np.diff(u)/2
    mid  = (u[:-1]+u[1:])/2
    xq   = np.r_[mid - half/np.sqrt(3), mid + half/np.sqrt(3)]
    wq   = np.r_[half, half]
    L    = np.sqrt(wq)[:,None]*BSpline(t, np.eye(nB), 3).derivative(2)(xq)
    Q, R = np.linalg.qr(np.r_[M@B, np.sqrt(lam)*L])
    Qd   = Q[:len(x)]
    # Omega = L^T L misses the straight lines, so its rank is nB-2.
    occam = -np.sum(np.log(np.abs(np.diag(R)))) + 0.5*(nB-2)*np.log(lam)
    return t, M, Qd, R, occam


def _logEvidence(My, QdTMy, occam):
    """
    Log evidence log p(y|lam) of a smoothing spline, up to a constant that does not depend on lam. The
    penalty is a Gaussian prior on the coefficients with precision lam Omega, flat in the straight-line
    directions, which Omega does not see. Integrating out c gives
        -2 log p(y|lam) = chi^2 + lam c^T Omega c + log det(X^T X) - log det+(lam Omega) + const,
    evaluated at the spline c, where det+ is the product of nonzero eigenvalues (REML). The first two terms
    are the minimum of |X c - b|^2, which is |M y|^2 - |Qd^T M y|^2, and log det(X^T X) = 2 sum log|R_ii|.
    """
    return -0.5*(My@My - np.sum(QdTMy**2,axis=-1)) + occam


def _lamGrid(xdata, edata, nlam):
    """
    Grid of nlam values of lam, uniform in log lam, for averaging smoothing splines. tr(H) falls
    monotonically with lam and does not depend on y, so we bracket the grid by the lam where tr(H) is
    close to the number of unique x (interpolation) and where it is close to 2 (straight line). This
    makes the grid independent of the units of x and y.
    """
    nu = len(np.unique(xdata))
    if nu < 3:
        logger.TBRaise('Averaging over lam needs at least 3 unique x.')

    def peffDiff(loglam,target):
        _, _, Qd, _, _ = _smoothingSystem(xdata, edata, np.exp(loglam))
        return np.sum(Qd**2) - target

    loglams = []
    for target in [nu-0.1, 2.01]:
        lo, hi = -5., 5.
        while peffDiff(lo,target) < 0:
            lo -= 10
        while peffDiff(hi,target) > 0:
            hi += 10
        loglams.append(brentq(peffDiff,lo,hi,args=(target,)))
    return np.exp(np.linspace(loglams[0],loglams[1],nlam))


def _averageWeights(y, M, Qd, occam):
    """
    Weights of the smoothing splines on the lam grid, proportional to their evidence, together with the
    log evidences. The arrays Qd and occam have the lam grid as their first axis.
    """
    My     = M @ y
    QdTMy  = np.einsum('knb,n->kb', Qd, My)
    logEv  = _logEvidence(My, QdTMy, occam)
    # The log evidence plays the role of -IC/2.
    return getModelWeights(-2*logEv), logEv, QdTMy


class TBSpline:

    """
    A fitted B-spline, as returned by getSpline, bootSpline, and posteriorSpline. It holds the data it was
    fit to, the tck of the fit in the format of splrep, and the kind of spline. Call it like a SciPy spline.
    For a smoothing spline, lam is the smoothing parameter, or 'average' for an average over lam, and peff
    is the effective number of parameters. A smoothing spline at fixed lam also knows its logEvidence.
    """

    def __init__(self,xdata,ydata,tck,edata=None,kind=None,natural=False,lam=None,peff=None,logEvidence=None):
        self.xspline     = np.copy(xdata)
        self.yspline     = np.copy(ydata)
        self.tck         = tck
        self.edata       = edata
        self.kind        = kind
        self.natural     = natural
        self.lam         = lam
        self.peff        = peff
        self.logEvidence = logEvidence
        # For a covariance matrix, weights holds 1/sqrt of its diagonal, for information.
        self.weights, self.cov = _parseEdata(edata)

    def get_nparams(self):
        """
        Number of free parameters. Natural boundary conditions remove 2. For a smoothing spline this
        is the effective number tr(H), which is not an integer. For a spline averaged over lam, it is
        the weighted average of tr(H).
        """
        if self.peff is not None:
            return self.peff
        nparams = len(self.get_knots()) - self.get_order() - 1
        if self.natural:
            nparams -= 2
        return nparams

    def get_AIC(self):
        """
        AIC of the fit, using the covariance matrix if edata is one. Requires edata.
        """
        if self.edata is None:
            logger.TBRaise('AIC needs edata.')
        cov = self.cov
        if cov is None:
            cov = np.diag(1/self.weights**2)
        return AIC(self.xspline, self.yspline, cov, self)

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
        return self.weights

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
    interpolation and averages over models.
    """

    hasAIC = True

    def fitter(self, xdata, edata=None):
        logger.TBRaise('Not implemented for',self)

    def fit(self, xdata, ydata, edata=None):
        """
        Fit this kind of spline to the data and return a TBSpline.
        """
        return self.fitter(xdata, edata).fit(np.asarray(ydata,dtype=float), warn=True)


class Interpolating(SplineKind):

    """
    Natural cubic interpolating spline, with a knot at every datum. It passes through the data, so edata
    are not used in the fit.
    """

    hasAIC = False

    def __repr__(self) -> str:
        return "Interpolating()"

    def fitter(self, xdata, edata=None):
        return _InterpolatingFitter(self, xdata, edata)


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

    def fitter(self, xdata, edata=None):
        return _RegressionFitter(self, xdata, edata)


class Smoothing(SplineKind):

    """
    Cubic smoothing spline. It minimizes (y - S)^T W (y - S) + lam int S''^2, with W as for Regression.
    For W = diag(w), this is the convention of scipy.interpolate.make_smoothing_spline. The minimizer is
    a natural cubic spline with a knot at every unique x, so there are no knots to choose. With
    lam='average', the result is instead an average of smoothing splines over nlam values of lam, uniform
    in log lam, weighted by their evidence. The grid runs from (almost) interpolating the data to (almost)
    a straight line. This needs edata.

    Args:
        lam (float or str, optional):
            Smoothing parameter > 0, or 'average'. Defaults to 'average'.
        nlam (int, optional):
            Number of values of lam for lam='average'. Defaults to 100.
    """

    def __init__(self, lam='average', nlam=100):
        # Set the attributes first, so that __repr__ works in error messages.
        self.lam      = lam
        self.nlam     = nlam
        self.averaged = isinstance(lam,str)
        if isinstance(lam,str):
            if lam != 'average':
                logger.TBRaise("lam must be a number or 'average'. lam =",lam)
            checkType('int',nlam=nlam)
            self.averaged = True
            self.hasAIC  = False
        else:
            checkType('real',lam=lam)
            if lam <= 0:
                logger.TBRaise('Need lam > 0. lam =',lam)
            self.averaged = False
            self.hasAIC  = True

    def __repr__(self) -> str:
        if self.averaged:
            return f"Smoothing(lam='average', nlam={self.nlam})"
        return f"Smoothing(lam={self.lam})"

    def fitter(self, xdata, edata=None):
        return _SmoothingFitter(self, xdata, edata)


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

    hasAIC = False

    def __init__(self, kmax=None, ngrid=None, natural=False, nupdates=10000, nburn=1000, every=10, nchains=1, 
                 seed=DEFAULTSEED):
        # Set the attributes first, so that __repr__ works in error messages.
        self.kmax     = kmax
        self.ngrid    = ngrid
        self.natural  = natural
        self.nupdates = nupdates
        self.nburn    = nburn
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
        checkType('int',nburn=nburn)
        checkType('int',every=every)
        checkType('int',nchains=nchains)
        checkType('int',seed=seed)
        if nburn >= nupdates:
            logger.TBRaise('Need nburn < nupdates. nburn, nupdates =',nburn,nupdates)

    def __repr__(self) -> str:
        return f"KnotAverage(kmax={self.kmax}, ngrid={self.ngrid}, natural={self.natural})"

    def fitter(self, xdata, edata=None):
        return _KnotAverageFitter(self, xdata, edata)


class _InterpolatingFitter:

    def __init__(self, kind, xdata, edata):
        self.kind  = kind
        self.x     = np.asarray(xdata,dtype=float)
        self.edata = edata

    def fit(self, y, warn=False):
        spl = make_interp_spline(self.x, y, k=3, bc_type='natural')
        # splrep pads the coefficients with order+1 zeros
        tck = (spl.t, np.r_[spl.c, np.zeros(4)], 3)
        return TBSpline(self.x, y, tck, edata=self.edata, kind=self.kind, natural=True)


def _lsqSystem(x, M, t, order, natural):
    """
    Set up the least-squares problem |M B N d - M y|^2 of a regression spline with full knot vector t by QR, 
    where c = N d. For a natural spline, the conditions S''(a) = S''(b) = 0 are C c = 0 and N is a basis of 
    the null space of C; otherwise N = 1. Returns N, Q, R, and whether the fit is underdetermined, i.e. 
    whether the Schoenberg-Whitney conditions fail.
    """
    nB = len(t) - order - 1
    B  = BSpline.design_matrix(x, t, order).toarray()
    if natural:
        C = BSpline(t, np.eye(nB), order).derivative(2)([t[0],t[-1]])
        N = null_space(C)
    else:
        N = np.eye(nB)
    X = M @ B @ N
    if X.shape[0] < X.shape[1]:
        return N, None, None, True
    Q, R = np.linalg.qr(X)
    d    = np.abs(np.diag(R))
    return N, Q, R, np.min(d) <= 1e-10*np.max(d)


class _RegressionFitter:

    def __init__(self, kind, xdata, edata):
        """
        Place the knots and set up the least-squares problem, see _lsqSystem.
        """
        self.kind  = kind
        self.x     = np.asarray(xdata,dtype=float)
        self.edata = edata
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
        for knot in kind.knots:
            knots.append(knot)
        knots = sorted(knots)
        if len(knots)>0 and knots[0]<xmin:
            logger.TBRaise("You can't put a knot to the left of the x-data. knots, min(xdata) = ",knots,xmin)
        if len(knots)>0 and knots[-1]>xmax:
            logger.TBRaise("You can't put a knot to the right of the x-data. knots, max(xdata) = ",knots,xmax)
        # Full knot vector: each endpoint repeated order+1 times around the interior knots, as in splrep.
        self.t = np.r_[[xmin]*(order+1), knots, [xmax]*(order+1)]
        self.M = _whitening(edata, len(self.x))
        self.N, self.Q, self.R, singular = _lsqSystem(self.x, self.M, self.t, order, kind.natural)
        if singular:
            logger.TBRaise('The fit is underdetermined. Each interval between knots needs enough data (Schoenberg-Whitney).')

    def fit(self, y, warn=False):
        d = solve_triangular(self.R, self.Q.T @ (self.M @ y))
        # splrep pads the coefficients with order+1 zeros
        tck = (self.t, np.r_[self.N @ d, np.zeros(self.kind.order+1)], self.kind.order)
        return TBSpline(self.x, y, tck, edata=self.edata, kind=self.kind, natural=self.kind.natural)


class _SmoothingFitter:

    def __init__(self, kind, xdata, edata):
        """
        Set up the system of _smoothingSystem for each lam: the given one, or the grid for lam='average'.
        Everything except the weights is linear in y.
        """
        self.kind  = kind
        self.x     = np.asarray(xdata,dtype=float)
        self.edata = edata
        if kind.averaged:
            if edata is None:
                logger.TBRaise("Averaging over lam needs edata, since the evidence depends on the errors.")
            self.lams = _lamGrid(self.x, edata, kind.nlam)
        else:
            self.lams = np.array([kind.lam])
        Qd, R, occam = [], [], []
        for lam in self.lams:
            t, M, Qdk, Rk, occamk = _smoothingSystem(self.x, edata, lam)
            Qd.append(Qdk)
            R.append(Rk)
            occam.append(occamk)
        self.t     = t
        self.M     = M
        self.Qd    = np.array(Qd)
        self.R     = np.array(R)
        self.occam = np.array(occam)
        # tr(H) = |Qd|^2, since H = Qd Qd^T
        self.peffs = np.sum(self.Qd**2,axis=(1,2))

    def _coeffs(self, y, warn):
        """
        Weights, log evidences, and coefficients at each lam.
        """
        w, logEv, QdTMy = _averageWeights(y, self.M, self.Qd, self.occam)
        if warn and self.kind.averaged and ((w[0] > 0.01) or (w[-1] > 0.01)):
            logger.warn('Weights at the ends of the lam grid are',w[0],'and',w[-1],'which may make the result',
                        'depend on the grid. Large weight at large lam means the data are consistent with a straight line.')
        c = []
        for k in range(len(self.lams)):
            c.append(solve_triangular(self.R[k], QdTMy[k]))
        return w, logEv, np.array(c)

    def fit(self, y, warn=False):
        """
        All splines share the knots, so their weighted average is again a spline, whose coefficients are the
        weighted average of the coefficients.
        """
        w, logEv, c = self._coeffs(y, warn)
        tck = (self.t, np.r_[w @ c, np.zeros(4)], 3)
        if self.kind.averaged:
            return TBSpline(self.x, y, tck, edata=self.edata, kind=self.kind, natural=True, lam='average',
                            peff=w @ self.peffs)
        return TBSpline(self.x, y, tck, edata=self.edata, kind=self.kind, natural=True, lam=self.lams[0],
                        peff=self.peffs[0], logEvidence=logEv[0])

    def details(self, y, xspl):
        """
        Everything about the average over lam for the data y: the weights, the smoothing spline at each lam,
        the averaged spline, and the weighted spreads over lam of the spline and of its xmax on xspl.
        """
        w, logEv, c = self._coeffs(y, False)
        Bspl = BSpline.design_matrix(xspl, self.t, 3).toarray()
        ylam, splinesLam = [], []
        for k in range(len(self.lams)):
            tck = (self.t, np.r_[c[k], np.zeros(4)], 3)
            splinesLam.append(TBSpline(self.x, y, tck, edata=self.edata, kind=Smoothing(lam=self.lams[k]),
                                       natural=True, lam=self.lams[k], peff=self.peffs[k], logEvidence=logEv[k]))
            ylam.append(Bspl @ c[k])
        ylam  = np.array(ylam)
        ybar  = w @ ylam
        xmaxs = xspl[np.argmax(ylam,axis=1)]
        det = {}
        det['weights']    = w
        det['logEv']      = logEv
        det['coeffs']     = c
        det['Bspl']       = Bspl
        det['splinesLam'] = splinesLam
        det['ybar']       = ybar
        det['ylamSpread'] = np.sqrt(w @ (ylam - ybar)**2)
        det['xmaxSpread'] = np.sqrt(w @ (xmaxs - w @ xmaxs)**2)
        return det


# The Metropolis-Hastings chain of KnotAverage. These functions use only arrays, loops, and numpy functions that numba
# supports, so that @compile turns them into machine code after numbaON(). The random numbers are drawn outside
# by TBRNG and passed in, so the chain is the same with and without numba.


@compile
def _bsplineDesign(x, t):
    """
    Dense design matrix B[i,j] = B_j(x[i]) of the cubic B-splines with the full knot vector t, whose 4 end knots
    coincide, by the Cox-de Boor recursion. Only the 4 B-splines that are nonzero on the interval of x[i] are
    computed.
    """
    n  = len(x)
    nB = len(t) - 4
    # Interval t[mu] <= x < t[mu+1]; the right end of the data belongs to the last interval.
    mu = np.searchsorted(t, x, side='right') - 1
    for i in range(n):
        if mu[i] > nB - 1:
            mu[i] = nB - 1
    # N[:,r] = B_{mu-d+r}(x) of degree d, built up from d=0
    N = np.zeros((n, 4))
    N[:,0] = 1.
    left  = np.zeros((n, 4))
    right = np.zeros((n, 4))
    for d in range(1, 4):
        left[:,d]  = x - t[mu+1-d]
        right[:,d] = t[mu+d] - x
        saved = np.zeros(n)
        for r in range(d):
            temp   = N[:,r]/(right[:,r+1] + left[:,d-r])
            N[:,r] = saved + right[:,r+1]*temp
            saved  = left[:,d-r]*temp
        N[:,d] = saved
    flat = np.zeros(n*nB)
    rows = np.arange(n)*nB
    for r in range(4):
        flat[rows + mu - 3 + r] = N[:,r]
    return flat.reshape((n, nB))


@compile
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


@compile
def _knotValid(idx, k, lt, le, nu):
    """
    Is there a unique datum strictly inside each interval of the configuration with grid indices idx[:k]?
    """
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


@compile
def _knotFit(x, M, My, grid, idx, k, natural):
    """
    Least-squares fit of the configuration with grid indices idx[:k] to the whitened data My. Returns whether
    the fit is determined, its AIC, and its B-spline coefficients.
    """
    t = np.zeros(k+8)
    t[:4]   = x.min()
    t[k+4:] = x.max()
    for i in range(k):
        t[4+i] = grid[idx[i]]
    B = _bsplineDesign(x, t)
    if natural:
        N = _naturalBasis(t)
    else:
        N = np.eye(k+4)
    X = np.ascontiguousarray(M @ B @ N)
    p = X.shape[1]
    if X.shape[0] < p:
        return False, 0., np.zeros(k+4)
    Q, R = np.linalg.qr(X)
    d = np.abs(np.diag(R))
    if d.min() <= 1e-10*d.max():
        return False, 0., np.zeros(k+4)
    QTMy = np.ascontiguousarray(Q.T) @ My
    chi2 = My @ My - QTMy @ QTMy
    c    = N @ np.linalg.solve(R, QTMy)
    return True, chi2 + 2*p, c


@compile
def _freePoint(idx, k, m):
    """
    The m-th grid index, counting from 0, that is not among the sorted indices idx[:k].
    """
    g = m
    for i in range(k):
        if idx[i] <= g:
            g += 1
        else:
            break
    return g


@compile
def _knotChain(x, M, My, grid, lt, le, nu, kmax, natural, idx0, k0, u, nburn, every, wmove, pbirth, pdeath):
    """
    Metropolis-Hastings chain over knot configurations, starting from idx0[:k0], with u[it,:] the 6 uniform
    random numbers of update it. Returns the grid indices (padded with -1), number of knots, AIC, and
    coefficients (padded with 0) of every saved configuration, and the proposal and acceptance counts of
    birth, death, and move.
    """
    G        = len(grid)
    nupdates = u.shape[0]
    nsave    = (nupdates - nburn - 1)//every + 1
    width    = max(kmax, 1)
    savedIdx = -np.ones((nsave, width), dtype=np.int64)
    savedK   = np.zeros(nsave, dtype=np.int64)
    savedAIC = np.zeros(nsave)
    savedC   = np.zeros((nsave, kmax+4))
    stats    = np.zeros((3, 2), dtype=np.int64)
    idx = np.zeros(width, dtype=np.int64)
    new = np.zeros(width, dtype=np.int64)
    k   = k0
    for i in range(k0):
        idx[i] = idx0[i]
    ok, aic, c = _knotFit(x, M, My, grid, idx, k, natural)
    isave = 0
    for it in range(nupdates):
        if k < kmax:
            b = pbirth
        else:
            b = 0.
        if k > 0:
            d = pdeath
        else:
            d = 0.
        propose  = False
        logRatio = 0.
        knew     = k
        if u[it,0] < b:
            move = 0
            # Birth at a random free grid point. Proposal ratio d_{k+1}/b_k; the grid factors cancel
            # against p(t|k), see the notes.
            g = _freePoint(idx, k, min(int(u[it,1]*(G-k)), G-k-1))
            j = 0
            while (j < k) and (idx[j] < g):
                new[j] = idx[j]
                j += 1
            new[j] = g
            for i in range(j, k):
                new[i+1] = idx[i]
            knew = k+1
            logRatio = np.log(pdeath/b)
            propose  = True
        elif u[it,0] < b + d:
            move = 1
            j = min(int(u[it,1]*k), k-1)
            for i in range(j):
                new[i] = idx[i]
            for i in range(j+1, k):
                new[i-1] = idx[i]
            knew = k-1
            logRatio = np.log(pbirth/d)
            propose  = True
        elif k > 0:
            move = 2
            j = min(int(u[it,1]*k), k-1)
            if u[it,2] < 0.5:
                # Relocate knot j to a random grid point not occupied by any knot.
                g = _freePoint(idx, k, min(int(u[it,3]*(G-k)), G-k-1))
                propose = True
            else:
                # Shift knot j by up to wmove grid points.
                step = 1 + min(int(u[it,3]*wmove), wmove-1)
                if u[it,4] < 0.5:
                    step = -step
                g = idx[j] + step
                propose = (g >= 0) and (g < G)
                for i in range(k):
                    if idx[i] == g:
                        propose = False
            if propose:
                # Remove knot j and insert g, keeping the indices sorted.
                m = 0
                inserted = False
                for i in range(k):
                    if i == j:
                        continue
                    if (not inserted) and (g < idx[i]):
                        new[m] = g
                        m += 1
                        inserted = True
                    new[m] = idx[i]
                    m += 1
                if not inserted:
                    new[m] = g
        else:
            move = 2
        stats[move,0] += 1
        if propose and _knotValid(new, knew, lt, le, nu):
            oknew, aicnew, cnew = _knotFit(x, M, My, grid, new, knew, natural)
            if oknew:
                logA = -0.5*(aicnew - aic) + logRatio
                if (logA >= 0) or ((logA > -700) and (u[it,5] < np.exp(logA))):
                    for i in range(knew):
                        idx[i] = new[i]
                    k   = knew
                    aic = aicnew
                    c   = cnew
                    stats[move,1] += 1
        if (it >= nburn) and ((it - nburn) % every == 0):
            for i in range(k):
                savedIdx[isave,i] = idx[i]
            savedK[isave]   = k
            savedAIC[isave] = aic
            savedC[isave,:k+4] = c
            isave += 1
    return savedIdx, savedK, savedAIC, savedC, stats


def _dataSeed(seed, ichain, y):
    """
    Seed for one chain from the seed of the kind, the chain index, and the bits of the data.
    """
    words = np.frombuffer(np.ascontiguousarray(y,dtype=np.float64).tobytes(),dtype=np.uint32)
    return [seed, ichain] + words.astype(np.int64).tolist()


class _KnotAverageFitter:

    # Probability to propose a birth (death) when k < kmax (k > 0); the rest are moves.
    PBIRTH = 1/3
    PDEATH = 1/3
    # Warn if the Geweke z of the AIC trace exceeds DRIFTZ, or if a chain saved fewer than MINSAVED configurations.
    DRIFTZ   = 3
    MINSAVED = 50

    def __init__(self, kind, xdata, edata):
        """
        Set up the grid of candidate knots, the largest number of knots, and the basis with a knot at every 
        grid point, in which the averaged spline is expressed.
        """
        if edata is None:
            logger.TBRaise('KnotAverage needs edata, since the AIC depends on the errors.')
        self.kind  = kind
        self.x     = np.asarray(xdata,dtype=float)
        self.edata = edata
        self.M     = _whitening(edata, len(self.x))
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
        self.tFull = np.r_[[self.xmin]*4, self.grid, [self.xmax]*4]
        self.z     = (self.tFull[1:-3] + self.tFull[2:-2] + self.tFull[3:-1])/3
        self.Pz    = np.linalg.inv(BSpline.design_matrix(self.z, self.tFull, 3).toarray())
        self._cacheKey = None
        self._cacheRun = None

    def _valid(self, idx):
        """
        Is there a unique datum strictly inside each interval of the configuration with grid indices idx?
        """
        if len(idx) == 0:
            return True
        if self.lt[idx[0]] < 2:
            return False
        if self.nu - self.le[idx[-1]] < 2:
            return False
        for i in range(len(idx)-1):
            if self.lt[idx[i+1]] - self.le[idx[i]] < 1:
                return False
        return True

    def _system(self, idx, My):
        """
        Fit of the configuration with grid indices idx to the whitened data My. Returns its AIC, knot vector, 
        and coefficients, or None if the configuration is not allowed.
        """
        if not self._valid(idx):
            return None
        t = np.r_[[self.xmin]*4, self.grid[idx], [self.xmax]*4]
        N, Q, R, singular = _lsqSystem(self.x, self.M, t, 3, self.kind.natural)
        if singular:
            return None
        QTMy = Q.T @ My
        chi2 = My @ My - QTMy @ QTMy
        c    = N @ solve_triangular(R, QTMy)
        return {'AIC':chi2 + 2*N.shape[1], 't':t, 'c':c, 'p':N.shape[1]}

    def _chain(self, My, rng, saved, stats, trace):
        """
        One Metropolis-Hastings chain for the whitened data My, run by _knotChain. Adds the saved configurations 
        to saved, a dict from the tuple of grid indices to [count, fit], the proposal and acceptance counts to 
        stats, and the AIC of every saved configuration, in order, to trace.
        """
        kind = self.kind
        # Start from up to 2 knots, evenly spaced in data index, moved to the nearest grid points.
        idx0 = np.zeros(max(self.kmax,1),dtype=np.int64)
        k0   = 0
        if min(2, self.kmax) > 0:
            start = np.unique(np.argmin(np.abs(np.subtract.outer(_even_knots(self.x,min(2,self.kmax)), self.grid)),axis=1))
            if self._system(start, My) is not None:
                k0 = len(start)
                idx0[:k0] = start
        # 6 uniform random numbers per update: move type, which knot or grid point, relocate or shift, where to,
        # direction of the shift, and acceptance.
        u = rng.random((kind.nupdates, 6))
        # Local moves shift a knot by up to wmove grid points.
        wmove = max(1, self.ngrid//20)
        savedIdx, savedK, savedAIC, savedC, counts = _knotChain(self.x, np.ascontiguousarray(self.M,dtype=float), 
            np.ascontiguousarray(My,dtype=float), self.grid, self.lt, self.le, self.nu, self.kmax, kind.natural, 
            idx0, k0, u, kind.nburn, kind.every, wmove, self.PBIRTH, self.PDEATH)
        imove = 0
        for move in ['birth','death','move']:
            stats[move][0] += counts[imove,0]
            stats[move][1] += counts[imove,1]
            imove += 1
        for isave in range(len(savedK)):
            k   = savedK[isave]
            key = tuple(savedIdx[isave,:k].tolist())
            trace.append(savedAIC[isave])
            if key in saved:
                saved[key][0] += 1
            else:
                t = np.r_[[self.xmin]*4, self.grid[list(key)], [self.xmax]*4]
                c = savedC[isave,:k+4].copy()
                if kind.natural:
                    p = k+2
                else:
                    p = k+4
                saved[key] = [1, {'AIC':savedAIC[isave], 't':t, 'c':c, 'p':p}]

    def _run(self, y):
        """
        Run the chains for the data y and pool their saved configurations. The result for the last y is 
        cached, since details needs it again.
        """
        key = np.asarray(y,dtype=float).tobytes()
        if key == self._cacheKey:
            return self._cacheRun
        My    = self.M @ y
        saved = {}
        stats  = {'birth':[0,0], 'death':[0,0], 'move':[0,0]}
        traces = []
        for ichain in range(self.kind.nchains):
            rng   = TBRNG(_dataSeed(self.kind.seed, ichain, y))
            trace = []
            self._chain(My, rng, saved, stats, trace)
            traces.append(np.array(trace))
        configs, counts, fits = [], [], []
        for cfg in saved:
            configs.append(self.grid[list(cfg)])
            counts.append(saved[cfg][0])
            fits.append(saved[cfg][1])
        counts = np.array(counts)
        run = {'configs':configs, 'counts':counts, 'fits':fits, 'stats':stats, 'traces':traces}
        self._cacheKey = key
        self._cacheRun = run
        return run

    def fit(self, y, warn=False):
        """
        The chain visits configurations with frequency proportional to their weight, so the average is the
        plain average over the saved configurations, weighted by how often each was saved.
        """
        run    = self._run(y)
        counts = run['counts']
        sz     = np.zeros(len(self.z))
        peff   = 0.
        kcount = np.zeros(self.kmax+1)
        for i in range(len(counts)):
            fit   = run['fits'][i]
            sz   += counts[i]*(BSpline.design_matrix(self.z, fit['t'], 3).toarray() @ fit['c'])
            peff += counts[i]*fit['p']
            kcount[len(run['configs'][i])] += counts[i]
        ntot  = np.sum(counts)
        sz   /= ntot
        peff /= ntot
        if warn and (kcount[self.kmax]/ntot > 0.01):
            logger.warn('The chain spends',kcount[self.kmax]/ntot,'of its time at kmax =',self.kmax,'knots. The data',
                        'may want more knots than allowed, or the configurations are fitting noise.')
        if warn:
            self._checkDrift(run['traces'])
        tck = (self.tFull, np.r_[self.Pz @ sz, np.zeros(4)], 3)
        return TBSpline(self.x, y, tck, edata=self.edata, kind=self.kind, natural=self.kind.natural, peff=peff)

    def _checkDrift(self, traces):
        """
        Geweke test of each chain: compare the mean AIC of the first 10% of its saved configurations with the
        mean AIC of the last 50%. The standard error of each mean comes from 5 batch means, which accounts for
        autocorrelation as long as a batch is longer than the autocorrelation time. A thermalized chain gives
        |z| of order 1; a chain still drifting after nburn gives a large |z|.
        """
        for ichain in range(len(traces)):
            trace = traces[ichain]
            if len(trace) < self.MINSAVED:
                logger.warn('Chain',ichain,'saved only',len(trace),'configurations, too few to check whether it',
                            'thermalized. Increase nupdates.')
                continue
            first = trace[:len(trace)//10]
            last  = trace[len(trace)//2:]
            se2   = _batchMeanError(first)**2 + _batchMeanError(last)**2
            diff  = np.mean(first) - np.mean(last)
            if se2 == 0:
                if diff == 0:
                    continue
                z = np.inf
            else:
                z = diff/np.sqrt(se2)
            if np.abs(z) > self.DRIFTZ:
                logger.warn('Chain',ichain,'may not be thermalized: the mean AIC of the first 10% of its saved',
                            'configurations differs from that of the last 50% by z =',round(z,1),'standard errors.',
                            'Increase nburn and nupdates.')

    def details(self, y, xspl):
        """
        The saved configurations for the data y, their counts, the histogram of k, the acceptance rates, and the
        weighted spreads over configurations of the spline and of its xmax on xspl.
        """
        run    = self._run(y)
        counts = run['counts']
        w      = counts/np.sum(counts)
        ycfg   = []
        kHist  = np.zeros(self.kmax+1,dtype=int)
        for i in range(len(counts)):
            fit = run['fits'][i]
            ycfg.append(BSpline.design_matrix(xspl, fit['t'], 3).toarray() @ fit['c'])
            kHist[len(run['configs'][i])] += counts[i]
        ycfg  = np.array(ycfg)
        ybar  = w @ ycfg
        xmaxs = xspl[np.argmax(ycfg,axis=1)]
        acceptance = {}
        for move in run['stats']:
            nprop, nacc = run['stats'][move]
            if nprop > 0:
                acceptance[move] = nacc/nprop
            else:
                acceptance[move] = 0.
        det = {}
        det['configs']     = run['configs']
        det['counts']      = counts
        det['kHistogram']  = kHist
        det['acceptance']  = acceptance
        det['ngrid']       = self.ngrid
        det['kmax']        = self.kmax
        det['yknotSpread'] = np.sqrt(w @ (ycfg - ybar)**2)
        det['xmaxSpread']  = np.sqrt(w @ (xmaxs - w @ xmaxs)**2)
        return det


def _batchMeanError(a, nbatch=5):
    """
    Error of the mean of the time series a from the spread of the means of nbatch consecutive batches.
    """
    means = []
    for batch in np.array_split(a, nbatch):
        means.append(np.mean(batch))
    return np.std(means, ddof=1)/np.sqrt(nbatch)


def _checkEdata(xdata, edata):
    if edata is None:
        return
    if np.ndim(edata) == 0:
        logger.TBRaise('edata must be a vector of errors or a covariance matrix. To choose a number of knots,',
                       'pass kind=Regression(nknots=...).')
    checkEqualLengths(xdata,edata)


def getSpline(xdata, ydata, edata=None, kind=None):
    """
    Fit a spline of the given kind to the data.

    Args:
        xdata (array-like)
        ydata (array-like)
        edata (array-like, optional):
            Error data, either a vector of errors or the covariance matrix of ydata. Defaults to None.
        kind (SplineKind, optional):
            Interpolating(), Regression(...), Smoothing(...), or KnotAverage(...). Defaults to Interpolating() without edata
            and Smoothing(), i.e. the average over lam, with edata.

    Returns:
        TBSpline
    """
    checkEqualLengths(xdata,ydata)
    _checkEdata(xdata,edata)
    if kind is None:
        if edata is None:
            kind = Interpolating()
        else:
            kind = Smoothing()
    return kind.fit(xdata, ydata, edata)


def getSplineErr(xdata, xspline, ydata, ydatae, kind=None):
    """
    Fit unweighted splines of the given kind to ydata-ydatae and ydata+ydatae, evaluate them at xspline,
    and return their midpoint and half-difference. This gives a smooth band, but it is not a statistical
    error propagation; for that use bootSpline. kind defaults to Interpolating().
    """
    if np.ndim(ydatae) != 1:
        logger.TBRaise('getSplineErr needs a vector of errors ydatae.')
    if kind is None:
        kind = Interpolating()
    spline_lower  = getSpline(xdata, ydata - ydatae, kind=kind)(xspline)
    spline_upper  = getSpline(xdata, ydata + ydatae, kind=kind)(xspline)
    spline_center = (spline_lower+spline_upper)/2
    spline_err    = (spline_upper-spline_lower)/2
    return spline_center, spline_err


def _medianBiasCorrect(f0, fBS):
    """
    Bias-corrected central value 2 f0 - median(fBS), with f0 the observable for the original data and fBS
    its bootstrap samples along axis 0. This is _biasCorrect of the bootstrap with the median in place of
    the mean.
    """
    return 2*np.asarray(f0) - std_median(fBS, axis=0)


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


def _kindSample(sample, fitter, xspl):
    """
    Fit one bootstrap sample of shape (1,ndata) for bootSpline. Return everything needed to rebuild the
    spline, plus the observables, as one flat array: yBS, coefficients, spline on xspl, the x on xspl where
    the spline is largest, and the number of parameters. This lives at module level so that it can be
    pickled for nproc>1.
    """
    yBS = sample[0]
    spl = fitter.fit(yBS)
    ys  = spl(xspl)
    return np.r_[yBS, spl.get_coeffs(), ys, xspl[np.argmax(ys)], spl.get_nparams()]


def bootSpline(xdata, ydata, edata, kind=None, numb_samples=300, nsupport=301, seed=DEFAULTSEED, nproc=1) -> dict:
    """
    Given xdata, ydata, edata, create a spline of the given kind. Use a Gaussian (parametric) bootstrap,
    drawing ydata from normal distributions of width edata, to propagate uncertainties of the data into an
    error band for the spline. Every bootstrap sample is fit with the same kind; e.g. for an average over
    lam, the weights are recomputed in every sample, so the error includes the uncertainty in lam. If edata
    is the covariance matrix of ydata, the samples are drawn from a multivariate normal, and the fits use
    the covariance matrix as well.

    Args:
        xdata (array-like)
        ydata (array-like)
        edata (array-like):
            Errors of ydata, or their covariance matrix.
        kind (SplineKind, optional):
            Kind of spline. Defaults to Smoothing(), i.e. the average over lam.
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
                the original data and f* from the samples; the error is dev_by_dist of the samples.
            xmax, xmaxe: the same for the x on xspl where the spline is largest
            xinflect, xinflecte, yinflect, yinflecte: the same for the inflection point with the largest |slope|.
                Samples without an inflection point are dropped, and all four are nan if splineMean has none.
            splineMean: the spline fit to the original data, not an average over samples
            splineBS: the spline of each bootstrap sample
            AIC: AIC of splineMean, unless the kind is Interpolating or averaged
        For Smoothing(lam='average') also:
            lams, peffs, weights, logEvidences: the grid of lam, the effective numbers of parameters, the
                weights, and the log evidences, for the original data.
            splinesLam: the smoothing spline at each lam for the original data.
            ysple_lam, xmaxe_lam: the weighted spreads over lam for the original data. They are diagnostics
                only: the bootstrap already contains the uncertainty in lam, so adding them to ysple or
                xmaxe would count it twice.
        For KnotAverage also, for the original data:
            knotConfigs, knotCounts: the distinct saved knot configurations and how often each was saved.
            kHistogram: how often each number of knots 0,...,kmax was saved.
            acceptance: acceptance rates of birth, death, and move proposals.
            ngrid, kmax: the number of candidate knots and the largest number of knots.
            ysple_knots, xmaxe_knots: the weighted spreads over configurations. Diagnostics only, like
                ysple_lam.
    """
    checkType("int",seed=seed)
    checkEqualLengths(xdata,ydata,edata)
    _checkEdata(xdata,edata)
    if kind is None:
        kind = Smoothing()
    y = np.asarray(ydata,dtype=float)
    weights, cov = _parseEdata(edata)
    std = 1/weights
    if cov is None:
        useCovariance = False
        Covariance    = None
    else:
        useCovariance = True
        Covariance    = [cov]
    xspl  = np.linspace(np.min(xdata),np.max(xdata),nsupport)
    ndata = len(y)

    fitter  = kind.fitter(xdata, edata)
    spline0 = fitter.fit(y, warn=True)
    t, _, k = spline0.tck
    ncoeff  = len(t)

    samples, _, _ = bootstr_from_gauss(_kindSample,data=[y],data_std_dev=[std],
                                       numb_samples=numb_samples,return_sample=True,seed=seed,
                                       args={'fitter':fitter,'xspl':xspl},nproc=nproc,
                                       useCovariance=useCovariance,Covariance=Covariance)
    samples = np.asarray(samples)

    # Rebuild the bootstrap splines without refitting.
    splines = []
    for sample in samples:
        yBS = sample[:ndata]
        c   = sample[ndata:ndata+ncoeff]
        if spline0.peff is None:
            peff = None
        else:
            peff = sample[-1]
        splines.append(TBSpline(xdata,yBS,(t,c,k),edata=edata,kind=kind,natural=spline0.natural,
                                lam=spline0.lam,peff=peff))

    res = {}
    res['splineBS']   = splines # in case you want spline functions at bootstrap level
    res['xspl']       = xspl
    # Central values are 2 f(original data) - median over samples, and errors are dev_by_dist over samples.
    # spline0 is the fit to the original data, which is what the bootstrap's zero-noise sample would give.
    ys0               = spline0(xspl)
    ysBS              = samples[:,ndata+ncoeff:-2]
    res['yspl']       = _medianBiasCorrect(ys0, ysBS)
    res['ysple']      = dev_by_dist(ysBS)
    res['splineMean'] = spline0
    # The second to last entry of each sample is that spline's xmax. Because of the bias correction, xmax can
    # lie between grid points.
    xmaxBS            = samples[:,-2]
    res['xmax']       = _medianBiasCorrect(xspl[np.argmax(ys0)], xmaxBS)
    res['xmaxe']      = dev_by_dist(xmaxBS)
    # The inflection point is found from the rebuilt splines, so samples without one can be dropped.
    xinf0, yinf0 = _steepestInflection(spline0, xspl)
    xinfBS, yinfBS = [], []
    for spline in splines:
        xinf, yinf = _steepestInflection(spline, xspl)
        if not np.isnan(xinf):
            xinfBS.append(xinf)
            yinfBS.append(yinf)
    if len(xinfBS) < len(splines):
        logger.warn(f'Dropped {len(splines)-len(xinfBS)} of {len(splines)} bootstrap samples without an inflection point.')
    if np.isnan(xinf0) or len(xinfBS) == 0:
        res['xinflect'], res['xinflecte'], res['yinflect'], res['yinflecte'] = np.nan, np.nan, np.nan, np.nan
    else:
        res['xinflect']  = _medianBiasCorrect(xinf0, xinfBS)
        res['xinflecte'] = dev_by_dist(xinfBS)
        res['yinflect']  = _medianBiasCorrect(yinf0, yinfBS)
        res['yinflecte'] = dev_by_dist(yinfBS)
    if kind.hasAIC:
        res['AIC'] = spline0.get_AIC()
    if isinstance(kind,KnotAverage):
        det = fitter.details(y, xspl)
        res['ysple_knots']  = det['yknotSpread']
        res['xmaxe_knots']  = det['xmaxSpread']
        res['knotConfigs']  = det['configs']
        res['knotCounts']   = det['counts']
        res['kHistogram']   = det['kHistogram']
        res['acceptance']   = det['acceptance']
        res['ngrid']        = det['ngrid']
        res['kmax']         = det['kmax']
    if isinstance(kind,Smoothing) and kind.averaged:
        det = fitter.details(y, xspl)
        res['ysple_lam']    = det['ylamSpread']
        res['xmaxe_lam']    = det['xmaxSpread']
        res['lams']         = fitter.lams
        res['peffs']        = fitter.peffs
        res['weights']      = det['weights']
        res['logEvidences'] = det['logEv']
        res['splinesLam']   = det['splinesLam']
    return res


def posteriorSpline(xdata, ydata, edata, kind=None, numb_samples=300, nsupport=301, seed=DEFAULTSEED) -> dict:
    """
    Smoothing spline and its Bayesian error band. The smoothing spline at fixed lam is the posterior mean
    of a Bayesian model, in which the curvature penalty is a Gaussian prior on the B-spline coefficients.
    Their posterior is Gaussian with covariance (B^T W B + lam Omega)^-1. With lam='average', we average
    over the lam grid, weighted by the evidence, i.e. integrate lam out with a flat prior in log lam. This
    is not a bootstrap: the data are used once, as measured, and the error band says how uncertain the true
    curve is, given the data and the assumption that it is smooth. Unlike the bootstrap error of bootSpline,
    this includes the uncertainty coming from the smoothing itself.

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
            lams, peffs, weights, logEvidences: the grid of lam, the effective numbers of parameters,
                the weights, and the log evidences.
            splinesLam: the smoothing spline at each lam, i.e. the posterior mean at fixed lam. Each
                has covPost, the posterior covariance of its coefficients.
    """
    checkType("int",seed=seed)
    checkType("int",numb_samples=numb_samples)
    checkEqualLengths(xdata,ydata,edata)
    _checkEdata(xdata,edata)
    if kind is None:
        kind = Smoothing()
    if not isinstance(kind,Smoothing):
        logger.TBRaise('posteriorSpline needs a Smoothing kind. Got',kind)
    y       = np.asarray(ydata,dtype=float)
    xspl    = np.linspace(np.min(xdata),np.max(xdata),nsupport)
    fitter  = kind.fitter(xdata, edata)
    spline0 = fitter.fit(y, warn=True)
    det     = fitter.details(y, xspl)
    w       = det['weights']
    Bspl    = det['Bspl']
    nlam    = len(fitter.lams)
    nB      = fitter.R.shape[2]

    # At fixed lam, the posterior covariance of the coefficients is A^-1 = R^-1 R^-T, so c_k + R^-1 z with
    # standard normal z is a draw from the posterior. Also accumulate the exact posterior variance on xspl.
    Rinv  = []
    ystat = np.zeros(len(xspl))
    for k in range(nlam):
        Rinvk = solve_triangular(fitter.R[k], np.eye(nB))
        Rinv.append(Rinvk)
        det['splinesLam'][k].covPost = Rinvk @ Rinvk.T
        ystat += w[k]*np.sum((Bspl @ Rinvk)**2,axis=1)
    ystat = np.sqrt(ystat)

    # Draw lam_k with probability w_k, then the coefficients
    rng = TBRNG(seed)
    ks  = rng.choice(nlam, size=numb_samples, p=w)
    z   = rng.standard_normal((numb_samples, nB))
    splines, xmaxs = [], []
    for i in range(numb_samples):
        k = ks[i]
        c = det['coeffs'][k] + Rinv[k] @ z[i]
        splines.append(TBSpline(xdata,y,(fitter.t,np.r_[c,np.zeros(4)],3),edata=edata,kind=kind,natural=True,
                                lam=spline0.lam,peff=spline0.peff))
        xmaxs.append(xspl[np.argmax(Bspl @ c)])
    xmaxs = np.array(xmaxs)

    res = {}
    res['xspl']          = xspl
    res['yspl']          = det['ybar']
    res['ysple_stat']    = ystat
    res['ysple_sys']     = det['ylamSpread']
    res['ysple']         = quadrature(np.array([ystat, det['ylamSpread']]),axis=0)
    res['splineMean']    = spline0
    res['splineSamples'] = splines
    res['xmax']          = std_median(xmaxs)
    res['xmaxe']         = dev_by_dist(xmaxs)
    res['xmaxe_sys']     = det['xmaxSpread']
    res['xmaxe_stat']    = np.sqrt(max(res['xmaxe']**2 - det['xmaxSpread']**2, 0))
    res['lams']          = fitter.lams
    res['peffs']         = fitter.peffs
    res['weights']       = w
    res['logEvidences']  = det['logEv']
    res['splinesLam']    = det['splinesLam']
    return res
