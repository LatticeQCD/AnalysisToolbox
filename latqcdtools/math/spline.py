#
# spline.py
#
# D. Clarke, Claude Code
#
# Splines with control over knot placement and endpoints, and their error bands. What kind of spline you
# get is set by a kind object: Interpolating, Regression, or KnotAverage. Pass it to
# getSpline for a single fit, or to bootSpline for a bootstrap error band.
#

import numpy as np
from scipy.interpolate import BSpline, make_interp_spline, splev
from scipy.linalg import null_space, solve_triangular
from scipy.optimize import brentq
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import AIC, whiteningMatrix, std_median, dev_by_dist
from latqcdtools.base.check import checkType, checkEqualLengths
from latqcdtools.base.initialize import TBRNG, DEFAULTSEED
from latqcdtools.base.speedify import compile
from latqcdtools.math.math import isMatrix, isVector, quadrature
from latqcdtools.statistics.bootstr import bootstr_from_gauss


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


class TBSpline:

    """
    A fitted B-spline, as returned by getSpline and bootSpline. It holds the data it was
    fit to, the tck of the fit in the format of splrep, and the kind of spline, which also says whether
    the spline is natural. Call it like a SciPy spline. For a spline averaged over knot configurations,
    peff is the average number of parameters.
    """

    def __init__(self,xdata,ydata,tck,kind,edata=None,peff=None):
        self.xspline     = np.copy(xdata)
        self.yspline     = np.copy(ydata)
        self.tck         = tck
        self.kind        = kind
        self.edata       = edata
        self.peff        = peff
        # For a covariance matrix, weights holds 1/sqrt of its diagonal, for information.
        self.weights, self.cov = _parseEdata(edata)

    @property
    def natural(self):
        return self.kind.natural

    def get_nparams(self):
        if self.peff is not None:
            return self.peff
        nparams = len(self.get_knots()) - self.get_order() - 1
        if self.natural:
            nparams -= 2
        return nparams

    def get_AIC(self):
        """
        AIC of the fit, using the covariance matrix if edata is one. Requires edata, and a kind with a
        meaningful AIC.
        """
        if not self.kind.hasAIC:
            logger.TBRaise('AIC is not meaningful for',self.kind)
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
    interpolation and averages over models. natural is True for kinds with zero curvature at the endpoints.
    """

    hasAIC  = True
    natural = False

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

    hasAIC  = False
    natural = True

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


class _Fitter:

    """
    Base class of the fitters. Every fitter has fit(y, warn), which returns a TBSpline; warn turns on
    warnings about the fit, where a fitter has any. Fitters of averaged kinds also have details(y, xspl),
    which returns the models of the average for the data y and their weighted spreads.
    """

    def __init__(self, kind, xdata, edata):
        self.kind  = kind
        self.x     = np.asarray(xdata,dtype=float)
        self.edata = edata

    def fit(self, y, warn=False):
        logger.TBRaise('Not implemented for',self.kind)


class _InterpolatingFitter(_Fitter):

    def fit(self, y, warn=False):
        spl = make_interp_spline(self.x, y, k=3, bc_type='natural')
        # splrep pads the coefficients with order+1 zeros
        tck = (spl.t, np.r_[spl.c, np.zeros(4)], 3)
        return TBSpline(self.x, y, tck, self.kind, edata=self.edata)


def _knotVector(xmin, xmax, interior, order=3):
    """
    Full knot vector as in splrep: each endpoint repeated order+1 times around the interior knots.
    """
    return np.r_[[xmin]*(order+1), interior, [xmax]*(order+1)]


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


class _RegressionFitter(_Fitter):

    def __init__(self, kind, xdata, edata):
        """
        Place the knots and set up the least-squares problem, see _lsqSystem.
        """
        super().__init__(kind, xdata, edata)
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
        self.t = _knotVector(xmin, xmax, knots, order)
        self.M = _whitening(edata, len(self.x))
        self.N, self.Q, self.R, singular = _lsqSystem(self.x, self.M, self.t, order, kind.natural)
        if singular:
            logger.TBRaise('The fit is underdetermined. Each interval between knots needs enough data (Schoenberg-Whitney).')

    def fit(self, y, warn=False):
        d = solve_triangular(self.R, self.Q.T @ (self.M @ y))
        # splrep pads the coefficients with order+1 zeros
        tck = (self.t, np.r_[self.N @ d, np.zeros(self.kind.order+1)], self.kind.order)
        return TBSpline(self.x, y, tck, self.kind, edata=self.edata)


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


def _nparams(k, natural):
    """
    Number of parameters of a cubic regression spline with k interior knots. Natural boundary conditions
    remove 2.
    """
    if natural:
        return k+2
    return k+4


def _dataSeed(seed, ichain, y):
    """
    Seed for one chain from the seed of the kind, the chain index, and the bits of the data.
    """
    words = np.frombuffer(np.ascontiguousarray(y,dtype=np.float64).tobytes(),dtype=np.uint32)
    return [seed, ichain] + words.astype(np.int64).tolist()


class _KnotAverageFitter(_Fitter):

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
        super().__init__(kind, xdata, edata)
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
        self.tFull = _knotVector(self.xmin, self.xmax, self.grid)
        # Clip, since rounding can put the end abscissae just outside [xmin, xmax], e.g. (3*1.6)/3 > 1.6.
        self.z     = np.clip((self.tFull[1:-3] + self.tFull[2:-2] + self.tFull[3:-1])/3, self.xmin, self.xmax)
        self.Pz    = np.linalg.inv(BSpline.design_matrix(self.z, self.tFull, 3).toarray())
        self._cacheKey = None
        self._cacheRun = None

    def _system(self, idx, My):
        """
        Fit of the configuration with grid indices idx to the whitened data My, by the same _knotValid and 
        _knotFit as the chain. Returns its AIC, knot vector, coefficients, and number of parameters, or None if 
        the configuration is not allowed.
        """
        idx = np.asarray(idx,dtype=np.int64)
        k   = len(idx)
        if not _knotValid(idx, k, self.lt, self.le, self.nu):
            return None
        ok, aic, c = _knotFit(self.x, np.ascontiguousarray(self.M,dtype=float), np.ascontiguousarray(My,dtype=float),
                              self.grid, idx, k, self.kind.natural)
        if not ok:
            return None
        t = _knotVector(self.xmin, self.xmax, self.grid[idx])
        return {'AIC':aic, 't':t, 'c':c, 'p':_nparams(k, self.kind.natural)}

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
                t = _knotVector(self.xmin, self.xmax, self.grid[list(key)])
                c = savedC[isave,:k+4].copy()
                saved[key] = [1, {'AIC':savedAIC[isave], 't':t, 'c':c, 'p':_nparams(k, kind.natural)}]

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
        # How often each number of knots 0,...,kmax was saved
        kHist  = np.zeros(self.kmax+1,dtype=int)
        for i in range(len(configs)):
            kHist[len(configs[i])] += counts[i]
        run = {'configs':configs, 'counts':counts, 'fits':fits, 'stats':stats, 'traces':traces, 'kHistogram':kHist}
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
        tck = (self.tFull, np.r_[self.Pz @ sz, np.zeros(4)], 3)
        return TBSpline(self.x, y, tck, self.kind, edata=self.edata, peff=peff)

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
        weighted spreads over configurations of the spline and its observables, see _modelSpreads.
        """
        run    = self._run(y)
        counts = run['counts']
        w      = counts/np.sum(counts)
        ycfg, splines = [], []
        for i in range(len(counts)):
            fit     = run['fits'][i]
            cfgKind = Regression(knots=run['configs'][i], natural=self.kind.natural)
            spline  = TBSpline(self.x, y, (fit['t'], np.r_[fit['c'], np.zeros(4)], 3), cfgKind, edata=self.edata)
            splines.append(spline)
            ycfg.append(spline(xspl))
        ycfg = np.array(ycfg)
        acceptance = {}
        for move in run['stats']:
            nprop, nacc = run['stats'][move]
            if nprop > 0:
                acceptance[move] = nacc/nprop
            else:
                acceptance[move] = 0.
        det = {}
        det['splines']     = splines
        det['weights']     = w
        det['configs']     = run['configs']
        det['counts']      = counts
        det['kHistogram']  = run['kHistogram']
        det['acceptance']  = acceptance
        det['ngrid']       = self.ngrid
        det['kmax']        = self.kmax
        det['spread']      = _modelSpreads(splines, ycfg, w, xspl)
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
            Interpolating(), Regression(...), or KnotAverage(...). Defaults to Interpolating().

    Returns:
        TBSpline
    """
    checkEqualLengths(xdata,ydata)
    _checkEdata(xdata,edata)
    if kind is None:
        kind = Interpolating()
    return kind.fit(xdata, ydata, edata)


def getSplineErr(xdata, xspline, ydata, ydatae, kind=None):
    """
    Fit unweighted splines of the given kind to ydata-ydatae and ydata+ydatae, evaluate them at xspline,
    and return their midpoint and half-difference. This gives a smooth band, but it is not a statistical
    error propagation; for that use bootSpline. kind defaults to Interpolating().
    """
    if np.ndim(ydatae) != 1:
        logger.TBRaise('getSplineErr needs a vector of errors ydatae.')
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


def _weightedSpread(w, a):
    """
    Spread sqrt(sum_k w_k (a_k - abar)^2) of a along axis 0, with abar = sum_k w_k a_k and normalized w.
    """
    return np.sqrt(w @ (a - w @ a)**2)


def _spreadOf(splines, w, f):
    """
    Weighted spread over the models splines, with weights w, of f(spline), which may be a number or an array.
    Models where f has a nan are dropped, and the remaining weights are renormalized. Returns nan if f has a
    nan for every model.
    """
    vals, wkeep = [], []
    for i in range(len(splines)):
        val = np.asarray(f(splines[i]),dtype=float)
        if not np.any(np.isnan(val)):
            vals.append(val)
            wkeep.append(w[i])
    if len(wkeep) == 0:
        return np.nan
    wkeep = np.array(wkeep)
    if np.sum(wkeep) < 0.99:
        logger.warn(f'Models where the quantity is undefined carry weight {1-np.sum(wkeep):.3f}; they are dropped from its systematic error.')
    return _weightedSpread(wkeep/np.sum(wkeep), np.array(vals))


def _maxOnGrid(ys, xspl):
    """
    For curves ys on xspl, one per row, the x on xspl where each is largest, and its value there.
    """
    imax = np.argmax(ys,axis=1)
    return xspl[imax], ys[np.arange(len(ys)),imax]


def _modelSpreads(splines, ys, w, xspl):
    """
    Systematic errors of an average over models (knot configurations): the weighted spreads over the
    models, with weights w, of each model's spline on xspl (ys holds them), its maximum (x and y on xspl), and
    its steepest inflection point. Models without an inflection point are dropped from that spread, see
    _spreadOf.
    """
    xmax, ymax = _maxOnGrid(ys, xspl)
    sp = {}
    sp['y']    = _weightedSpread(w, ys)
    sp['xmax'] = _weightedSpread(w, xmax)
    sp['ymax'] = _weightedSpread(w, ymax)
    def inflection(spline):
        return _steepestInflection(spline, xspl)
    spread = _spreadOf(splines, w, inflection)
    if np.all(np.isnan(spread)):
        sp['xinflect'], sp['yinflect'] = np.nan, np.nan
    else:
        sp['xinflect'], sp['yinflect'] = spread[0], spread[1]
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
    return _spreadOf(res['splinesModel'], res['weightsModel'], f)


def _kindSample(sample, fitter):
    """
    Fit one bootstrap sample of shape (1,ndata) for bootSpline. The bootstrap needs numbers, so pack what
    is needed to rebuild the spline into one flat array: yBS, the coefficients, and the number of parameters.
    _unpackSample undoes this. This lives at module level so that it can be pickled for nproc>1.
    """
    yBS = sample[0]
    spl = fitter.fit(yBS)
    return np.r_[yBS, spl.get_coeffs(), spl.get_nparams()]


def _unpackSample(sample, ndata):
    """
    yBS, the coefficients, and the number of parameters from an array packed by _kindSample.
    """
    return sample[:ndata], sample[ndata:-1], sample[-1]


def _setErrors(res, name, ename, f0, fBS, sys):
    """
    Put the bias-corrected central value of an observable into res[name], and its statistical error over
    the bootstrap samples fBS, its systematic error sys, and their sum in quadrature into res[ename+'_stat'],
    res[ename+'_sys'], and res[ename].
    """
    stat = dev_by_dist(fBS)
    res[name]          = _medianBiasCorrect(f0, fBS)
    res[ename+'_stat'] = stat
    res[ename+'_sys']  = sys
    res[ename]         = quadrature(np.array([stat, sys]),axis=0)


def bootSpline(xdata, ydata, edata, kind=None, numb_samples=300, nsupport=301, seed=DEFAULTSEED, nproc=1) -> dict:
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
    checkType("int",seed=seed)
    checkEqualLengths(xdata,ydata)
    if edata is None:
        logger.TBRaise('bootSpline draws its samples from edata, so it needs edata. For a spline without errors, use getSpline.')
    _checkEdata(xdata,edata)
    if kind is None:
        kind = Interpolating()
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

    samples, _, _ = bootstr_from_gauss(_kindSample,data=[y],data_std_dev=[std],
                                       numb_samples=numb_samples,return_sample=True,seed=seed,
                                       args={'fitter':fitter},nproc=nproc,
                                       useCovariance=useCovariance,Covariance=Covariance)
    samples = np.asarray(samples)

    # Rebuild the bootstrap splines without refitting.
    splines = []
    for sample in samples:
        yBS, c, nparams = _unpackSample(sample, ndata)
        if spline0.peff is None:
            peff = None
        else:
            peff = nparams
        splines.append(TBSpline(xdata,yBS,(t,c,k),kind,edata=edata,peff=peff))

    # Systematic errors from the spread over models for the original data. Kinds that are not averages have none.
    if isinstance(kind,KnotAverage):
        det    = fitter.details(y, xspl)
        spread = det['spread']
    else:
        det    = None
        spread = {'y':np.zeros(len(xspl)), 'xmax':0., 'ymax':0., 'xinflect':0., 'yinflect':0.}

    res = {}
    res['splineBS']   = splines # in case you want spline functions at bootstrap level
    res['xspl']       = xspl
    res['splineMean'] = spline0
    # Central values are 2 f(original data) - median over samples. The statistical errors are dev_by_dist over 
    # samples, and the totals add the systematic ones in quadrature. spline0 is the fit to the original data,
    # which is what the bootstrap's zero-noise sample would give.
    ys0     = spline0(xspl)
    xmax0, ymax0 = _maxOnGrid(np.array([ys0]), xspl)
    ysBS    = []
    for spline in splines:
        ysBS.append(spline(xspl))
    ysBS    = np.array(ysBS)
    xmaxBS, ymaxBS = _maxOnGrid(ysBS, xspl)
    # Because of the bias correction, xmax can lie between grid points.
    _setErrors(res, 'yspl', 'ysple', ys0, ysBS, spread['y'])
    _setErrors(res, 'xmax', 'xmaxe', xmax0[0], xmaxBS, spread['xmax'])
    _setErrors(res, 'ymax', 'ymaxe', ymax0[0], ymaxBS, spread['ymax'])
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
        for key in ['xinflect', 'xinflecte', 'xinflecte_stat', 'xinflecte_sys',
                    'yinflect', 'yinflecte', 'yinflecte_stat', 'yinflecte_sys']:
            res[key] = np.nan
    else:
        _setErrors(res, 'xinflect', 'xinflecte', xinf0, xinfBS, spread['xinflect'])
        _setErrors(res, 'yinflect', 'yinflecte', yinf0, yinfBS, spread['yinflect'])
    if kind.hasAIC:
        res['AIC'] = spline0.get_AIC()
    if isinstance(kind,KnotAverage):
        res['splinesModel'] = det['splines']
        res['weightsModel'] = det['weights']
        res['knotConfigs']  = det['configs']
        res['knotCounts']   = det['counts']
        res['kHistogram']   = det['kHistogram']
        res['acceptance']   = det['acceptance']
        res['ngrid']        = det['ngrid']
        res['kmax']         = det['kmax']
    return res
