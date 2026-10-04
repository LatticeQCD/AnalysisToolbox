# 
# spline.py                                                               
# 
# D. Clarke, Claude Code
# 
# Wrappers for SciPy splines that give control over knot placement and endpoints. Depending on the arguments,
# getSpline gives an interpolating natural cubic spline, a least-squares (regression) spline with fixed knots,
# or a smoothing spline with a curvature penalty. 
# 

import numpy as np
from scipy.interpolate import CubicSpline, BSpline, splrep, splev
from scipy.linalg import null_space, solve_triangular, cholesky
from scipy.optimize import brentq
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import AICc
from latqcdtools.base.check import checkType, checkEqualLengths
from latqcdtools.base.initialize import TBRNG, DEFAULTSEED
from latqcdtools.math.math import isMatrix, isVector
from latqcdtools.statistics.statistics import std_median, dev_by_dist
from latqcdtools.statistics.bootstr import bootstr_from_gauss
from latqcdtools.math.num_deriv import diff_deriv


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
        logger.TBRaise('Expected edata with ndim < 3. Got ndim =',edata.ndim)


def _smoothingSystem(xdata, edata, lam):
    """
    Set up and solve by QR the least-squares problem of a cubic smoothing spline, which minimizes 
    (y - S)^T W (y - S) + lam int S''^2 with W = diag(1/edata^2), or W = edata^-1 if edata is a covariance
    matrix. The minimizer over all twice-differentiable functions is a natural cubic spline with a knot at 
    every unique x, so we work in that B-spline basis. With S = B c and int S''^2 = |L c|^2, the problem 
    is |X c - b|^2 with X = (M B, sqrt(lam) L) and b = (M y, 0) stacked, where M^T M = W. For diagonal W, 
    M = sqrt(W), and otherwise M = C^-1 with C the Cholesky factor of the covariance. QR, X = Q R, stays 
    accurate for large lam, unlike the normal equations. None of this depends on y. 

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
    weights, cov = _parseEdata(edata)
    if cov is not None:
        M = solve_triangular(cholesky(cov, lower=True), np.eye(len(x)), lower=True)
    elif weights is None:
        M = np.eye(len(x))
    else:
        M = np.diag(weights)
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


class TBSpline:

    """
    Least-squares (regression) spline with fixed interior knots. Without natural, prepares a splrep
    with task=-1 and wraps it with splev. With natural, does the least-squares fit itself in the
    subspace of cubic B-splines with zero curvature at the endpoints. With lam, it is instead a cubic
    smoothing spline with a knot at every datum, and knots and natural are ignored. edata may be a 
    vector of errors or a covariance matrix. Only the smoothing spline can take a covariance matrix.
    """

    def __init__(self,xdata,ydata,edata=None,knots=None,order=3,natural=False,lam=None):

        self.xspline     = np.copy(xdata)
        self.yspline     = np.copy(ydata)
        self.natural     = natural
        self.lam         = lam
        self.peff        = None
        self.logEvidence = None
        # For a covariance matrix, weights holds 1/sqrt of its diagonal, for information.
        self.weights, self.cov = _parseEdata(edata)
        self.edata       = edata

        if (self.cov is not None) and (lam is None):
            logger.TBRaise('Only smoothing splines can take a covariance matrix.')
        if self.weights is None:
            smooth = None
        else:
            # Ignored by splrep when task=-1.
            smooth = len(self.weights)

        if lam is not None:
            # The minimizer is natural automatically.
            self.natural = True
            self.tck, self.peff, self.logEvidence = self._smoothingFit(lam)
        elif natural:
            self.tck = self._naturalFit(knots,order)
        else:
            self.tck = splrep(self.xspline, self.yspline, t=knots, k=order, w=self.weights, s=smooth, task=-1)

    @classmethod
    def fromTck(cls,xdata,ydata,tck,edata=None,natural=False,lam=None,peff=None):
        """
        Build a TBSpline from known tck without fitting, e.g. from coefficients found elsewhere. For a
        smoothing spline, also pass lam and peff, which do not depend on ydata.
        """
        spline = cls.__new__(cls)
        spline.xspline     = np.copy(xdata)
        spline.yspline     = np.copy(ydata)
        spline.natural     = natural
        spline.lam         = lam
        spline.peff        = peff
        spline.logEvidence = None
        spline.weights, spline.cov = _parseEdata(edata)
        spline.edata       = edata
        spline.tck         = tck
        return spline

    def _naturalFit(self,knots,order):
        """
        The conditions S''(a) = S''(b) = 0 are linear in the B-spline coefficients c, i.e. C c = 0.
        Write c = N d with N a basis of the null space of C, then do a weighted linear least-squares
        fit for d. Returns tck in the same format as splrep.
        """
        x = np.asarray(self.xspline,dtype=float)
        y = np.asarray(self.yspline,dtype=float)
        if knots is None:
            knots = []
        # Full knot vector: each endpoint repeated order+1 times around the interior knots, as in splrep.
        t = np.r_[[x[0]]*(order+1), knots, [x[-1]]*(order+1)]
        # Number of B-spline basis functions, i.e. unconstrained parameters.
        nB = len(t) - order - 1
        # Design matrix B[i,j] = B_j(x_i), so the spline at the data is B c.
        B  = BSpline.design_matrix(x, t, order).toarray()
        # Constraint matrix C[:,j] = (B_j''(a), B_j''(b)), so S''(a) = S''(b) = 0 reads C c = 0.
        C  = np.zeros((2,nB))
        for j in range(nB):
            Bj = BSpline(t, np.eye(nB)[j], order)
            C[:,j] = Bj.derivative(2)([x[0],x[-1]])
        # Orthonormal basis of the null space of C; any c = N d satisfies the constraints.
        N  = null_space(C)
        if self.weights is None:
            w = np.ones(len(x))
        else:
            w = self.weights
        # Weighted least squares for d: scale rows of B N and y by w = 1/sigma.
        d  = np.linalg.lstsq((B*w[:,None])@N, y*w, rcond=None)[0]
        # splrep pads the coefficients with order+1 zeros
        return (t, np.r_[N@d, np.zeros(order+1)], order)

    def _smoothingFit(self,lam):
        """
        Cubic smoothing spline, see _smoothingSystem. Returns tck in the same format as splrep, the
        effective number of parameters tr(H), where H is the hat matrix, and the log evidence (up to a 
        constant that does not depend on lam). Since H = Qd Qd^T, tr(H) = |Qd|^2.
        """
        y = np.asarray(self.yspline,dtype=float)
        t, M, Qd, R, occam = _smoothingSystem(self.xspline, self.edata, lam)
        My    = M @ y
        QdTMy = Qd.T @ My
        c     = solve_triangular(R, QdTMy)
        peff  = np.sum(Qd**2)
        # splrep pads the coefficients with order+1 zeros
        return (t, np.r_[c, np.zeros(4)], 3), peff, _logEvidence(My, QdTMy, occam)

    def get_nparams(self):
        """
        Number of free parameters. Natural boundary conditions remove 2. For a smoothing spline this
        is the effective number tr(H), which is not an integer. For a spline averaged over lam, it is
        the weighted average of tr(H).
        """
        if self.lam is not None:
            return self.peff
        nparams = len(self.get_knots()) - self.get_order() - 1
        if self.natural:
            nparams -= 2
        return nparams

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



def getSpline(xdata, ydata, num_knots=None, edata=None, order=3, rand=False, fixedKnots=None, 
              getAICc=False, natural=False, seed=None, lam=None):
    """ 
    This is a wrapper that calls SciPy spline methods, depending on your needs. If natural=True and
    edata=None, it returns a natural cubic interpolating spline from scipy.interpolate.CubicSpline,
    with a knot at every data point. Otherwise it returns a TBSpline, i.e. a (weighted, if edata are
    given) least-squares B-spline fit with fixed knots. If natural=True and edata are provided, this
    least-squares spline has zero curvature at the endpoints. If lam is given, it returns a cubic 
    smoothing spline instead, which has a knot at every unique x and is natural automatically. 

    Args:
        xdata (array-like)
        ydata (array-like)
        num_knots (int):
            The number of interior knots, including fixedKnots. Not used for the interpolating spline.
        edata (array-like, optional): 
            Error data, either a vector of errors or the covariance matrix of ydata. A covariance matrix
            requires lam. Defaults to None.
        order (int, optional):
            Degree of the spline polynomials (SciPy's k). Defaults to 3.
        rand (bool, optional): 
            Use randomly placed knots? Defaults to False.
        fixedKnots (list, optional):
            List of user-specified knots. These count toward num_knots. Defaults to None.
        getAICc (bool, optional): 
            Return corrected Akaike information criterion? Requires edata. Defaults to False.
        natural (bool, optional): 
            Try a natural (zero curvature at the endpoints) cubic spline. Defaults to False. 
        seed (int, optional):
            Seed for the random knots when rand=True. Defaults to None.
        lam (float, optional):
            Smoothing parameter of a smoothing spline, which minimizes sum_i w_i (y_i - S(x_i))^2 + lam int S''^2
            with w_i = 1/edata_i^2, or 1 without edata. This is the convention of 
            scipy.interpolate.make_smoothing_spline. If edata is a covariance matrix, the first term becomes
            (y - S)^T edata^-1 (y - S). The log evidence of the fit is stored in the logEvidence attribute.
            Defaults to None, i.e. no smoothing spline.

    Returns:
        callable spline object
        AICc (optionally)
    """

    if len(xdata) != len(ydata):
        logger.TBRaise('len(xdata), len(ydata) =',len(xdata),len(ydata))
    if natural and (order != 3):
        logger.TBRaise("Natural splines have order=3 by definition.")
    if getAICc and (edata is None):
        logger.TBRaise("getAICc requires edata.")

    if lam is not None:
        if (num_knots is not None) or (fixedKnots is not None) or rand:
            logger.TBRaise('A smoothing spline has a knot at every unique x. Do not pass num_knots, fixedKnots, or rand.')
        if order != 3:
            logger.TBRaise('Smoothing splines are only implemented for order=3.')
        if lam <= 0:
            logger.TBRaise('Need lam > 0. lam =',lam)
        spline = TBSpline(xdata, ydata, edata=edata, order=order, lam=lam)

    elif natural and (edata is None): 
        if num_knots is not None:
            logger.TBRaise('Natural spline without errors interpolates, with a knot at every data point. Do not pass num_knots.')
        spline = CubicSpline(x=xdata,y=ydata,bc_type='natural')

    else:
        checkType('int',num_knots=num_knots)
        nknots = num_knots
        if fixedKnots is not None:
            if type(fixedKnots) is not list:
                logger.TBRaise("knots must be specified as a list.")
            nknots -= len(fixedKnots)
            if nknots < 0:
                logger.TBRaise("len(fixedKnots)",len(fixedKnots),"exceeds num_knots",num_knots)
        if nknots>0:
            if rand:
                knots = _random_knots(xdata,nknots,SEED=seed)
            else:
                knots = _even_knots(xdata,nknots)
        else:
            knots=[]
        if fixedKnots is not None:
            for knot in fixedKnots:
                knots.append(knot)
        knots = sorted(knots)
        if len(knots)>0 and knots[0]<xdata[0]:
            logger.TBRaise("You can't put a knot to the left of the x-data. knots, xdata[0] = ",knots,xdata[0])
        if len(knots)>0 and knots[-1]>xdata[-1]:
            logger.TBRaise("You can't put a knot to the right of the x-data. knots, xdata[-1] = ",knots,xdata[-1])
        spline = TBSpline(xdata, ydata, edata=edata, knots=knots, order=order, natural=natural)

    if getAICc:
        _, cov = _parseEdata(edata)
        if cov is None:
            cov = np.diag(np.asarray(edata)**2)
        return spline, AICc(xdata, ydata, cov, spline)
    else:
        return spline


def getSplineErr(xdata, xspline, ydata, ydatae, num_knots=None, order=3, rand=False, fixedKnots=None, natural=False,
                 seed=None):
    """ 
    Fit unweighted least-squares splines to ydata-ydatae and ydata+ydatae, evaluate them at xspline,
    and return their midpoint and half-difference. This gives a smooth band, but it is not a
    statistical error propagation; for that use bootSpline. 
    """
    if np.ndim(ydatae) != 1:
        logger.TBRaise('getSplineErr needs a vector of errors ydatae.')
    # Both splines must use the same random knots.
    if rand and (seed is None):
        seed = int(TBRNG().integers(2**32))
    spline_lower  = getSpline(xdata, ydata - ydatae, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural, seed=seed)(xspline)
    spline_upper  = getSpline(xdata, ydata + ydatae, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural, seed=seed)(xspline)
    spline_center = (spline_lower+spline_upper)/2 
    spline_err    = (spline_upper-spline_lower)/2 
    return spline_center, spline_err


def _splineSample(sample, xdata, edata, xspl, num_knots, order, fixedKnots, natural, lam):
    """
    Fit one bootstrap sample of shape (1,ndata) for bootSpline. Return everything needed to rebuild the
    spline, plus the observables, as one flat array: yBS, coefficients, spline on xspl, and the x on
    xspl where the spline is largest. This lives at module level so that it can be pickled for nproc>1.
    """
    yBS = sample[0]
    spl = getSpline(xdata=xdata,ydata=yBS,num_knots=num_knots,edata=edata,order=order,
                    fixedKnots=fixedKnots,natural=natural,lam=lam)
    ys  = spl(xspl)
    return np.r_[yBS, spl.get_coeffs(), ys, xspl[np.argmax(ys)]]


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
    # Weights below exp(-70) ~ 1e-30 relative to the largest are negligible. Setting them to zero avoids
    # underflow, which the toolbox treats as an error.
    dLogEv = logEv - np.max(logEv)
    w      = np.zeros(len(logEv))
    keep   = dLogEv > -70
    w[keep] = np.exp(dLogEv[keep])
    return w/np.sum(w), logEv, QdTMy


def _averageSample(sample, M, Qd, R, occam, Bspl, xspl):
    """
    Average the smoothing splines on the lam grid for one bootstrap sample of shape (1,ndata), recomputing
    the weights. All splines share the knots, so the average is again a spline, whose coefficients are 
    the weighted average of the coefficients. Return as one flat array: yBS, averaged coefficients, the 
    averaged spline on xspl, and the x on xspl where it is largest. Module level, so it can be pickled.
    """
    yBS = sample[0]
    w, _, QdTMy = _averageWeights(yBS, M, Qd, occam)
    cbar = np.zeros(R.shape[1])
    for k in range(len(w)):
        cbar += w[k]*solve_triangular(R[k], QdTMy[k])
    ys = Bspl @ cbar
    return np.r_[yBS, cbar, np.zeros(4), ys, xspl[np.argmax(ys)]]


def bootSpline(xdata, ydata, edata, num_knots=None, order=3, fixedKnots=None, natural=False, numb_samples=300, 
               nsupport=301, seed=DEFAULTSEED, nproc=1, lam=None, nlam=100) -> dict:
    """
    Given xdata, ydata, edata, create a spline. Use a Gaussian (parametric) bootstrap, drawing ydata
    from normal distributions of width edata, to propagate uncertainties of the data into an error band
    for the spline. Gives back a dictionary whose xspl, yspl, and ysple entries can be used to plot a
    spline with error bars. yspl and xmax are bias-corrected central values and ysple and xmaxe their
    errors, as in bootstr_from_gauss. splineMean is the spline fit to the original data, not an average,
    and AICc is the AICc of that fit. The knots are the same for every bootstrap sample. nproc is
    passed to bootstr_from_gauss. With lam, the spline is a smoothing spline, see getSpline. If edata 
    is the covariance matrix of ydata, the samples are drawn from a multivariate normal. For a smoothing
    spline, the fits and the AICc then use the covariance matrix as well. Otherwise they are weighted by 
    the square root of its diagonal.

    With lam='average', the result is instead an average of smoothing splines over nlam values of lam, 
    weighted by their evidence. The lam grid is uniform in log lam and runs from (almost) interpolating 
    the data to (almost) a straight line. In every bootstrap sample the weights are recomputed, so ysple 
    and xmaxe include the uncertainty in lam. splineMean is the averaged spline for the original data. 
    You also get lams, the effective numbers of parameters peffs, the weights, and the logEvidences on 
    the grid, and splinesLam, the smoothing spline at each lam for the original data. ysple_lam and 
    xmaxe_lam are the weighted spreads over lam for the original data. They are diagnostics only: the 
    bootstrap already contains the uncertainty in lam, so adding them to ysple would count it twice. 
    For an error band that includes the uncertainty from the smoothing itself, see posteriorSpline.
    """
    checkType("int",seed=seed)
    checkEqualLengths(xdata,ydata,edata)
    edata = np.asarray(edata,dtype=float)
    weights, cov = _parseEdata(edata)
    std = 1/weights
    if cov is None:
        useCovariance = False
        Covariance    = None
    else:
        useCovariance = True
        Covariance    = [cov]
    # Weights of the fits
    if (lam is None) and (cov is not None):
        fitEdata = std
    else:
        fitEdata = edata
    xspl  = np.linspace(np.min(xdata),np.max(xdata),nsupport)
    ndata = len(ydata)

    if isinstance(lam,str):
        if lam != 'average':
            logger.TBRaise("lam must be a number or 'average'. lam =",lam)
        return _bootSplineAverage(xdata, ydata, fitEdata, std, useCovariance, Covariance, xspl, 
                                  num_knots, fixedKnots, order, numb_samples, seed, nproc, nlam)

    # The AICc must come from the original data. A bootstrap sample scatters around ydata, which already
    # scatters around the truth, so its chi^2 is about twice as large. 
    spline0, AICc0 = getSpline(xdata=xdata,ydata=ydata,edata=fitEdata,num_knots=num_knots,order=order,
                               fixedKnots=fixedKnots,natural=natural,getAICc=True,lam=lam)
    t, _, k = spline0.tck
    ncoeff  = len(t)

    fitArgs = {'xdata':xdata,'edata':fitEdata,'xspl':xspl,'num_knots':num_knots,'order':order,
               'fixedKnots':fixedKnots,'natural':natural,'lam':lam}
    samples, mean, err = bootstr_from_gauss(_splineSample,data=[ydata],data_std_dev=[std],
                                            numb_samples=numb_samples,return_sample=True,seed=seed,
                                            args=fitArgs,nproc=nproc,useCovariance=useCovariance,
                                            Covariance=Covariance)

    # Rebuild the bootstrap splines without refitting.
    splines = [] 
    for sample in samples:
        yBS = sample[:ndata]
        c   = sample[ndata:ndata+ncoeff]
        splines.append(TBSpline.fromTck(xdata,yBS,(t,c,k),edata=fitEdata,natural=spline0.natural,
                                        lam=lam,peff=spline0.peff))

    res = {}
    res['AICc']       = AICc0   # inf if nparams >= ndata-1
    res['splineBS']   = splines # in case you want spline functions at bootstrap level
    res['xspl']       = xspl
    res['yspl']       = mean[ndata+ncoeff:-1]
    res['ysple']      = err[ndata+ncoeff:-1]
    res['splineMean'] = spline0
    # The last entry of each sample is that spline's xmax, so mean[-1] and err[-1] are the bias-corrected
    # central value and error of xmax over the samples. Hence xmax can lie between grid points.
    res['xmax']       = mean[-1]
    res['xmaxe']      = err[-1]
                        
    return res


def _lamAverage(xdata, y, edata, xspl, nlam):
    """
    Everything about the average of smoothing splines over lam that concerns the original data y. Returns
    a dict with the grid, the systems of _smoothingSystem at each lam, the weights, the smoothing spline 
    at each lam, the averaged spline, and the weighted spreads over lam of the spline and of its xmax.
    """
    lams = _lamGrid(xdata, edata, nlam)

    # Everything except the weights is linear in y, so set up each lam once.
    Qd, R, occam = [], [], []
    for lam in lams:
        t, M, Qdk, Rk, occamk = _smoothingSystem(xdata, edata, lam)
        Qd.append(Qdk)
        R.append(Rk)
        occam.append(occamk)
    Qd, R, occam = np.array(Qd), np.array(R), np.array(occam)
    peffs = np.sum(Qd**2,axis=(1,2))
    Bspl  = BSpline.design_matrix(xspl, t, 3).toarray()

    w, logEv, QdTMy = _averageWeights(y, M, Qd, occam)
    if (w[0] > 0.01) or (w[-1] > 0.01):
        logger.warn('Weights at the ends of the lam grid are',w[0],'and',w[-1],'which may make the result',
                    'depend on the grid. Large weight at large lam means the data are consistent with a straight line.')

    # Smoothing spline at each lam, and their average
    ylam, splinesLam = [], []
    cbar = np.zeros(R.shape[2])
    for k in range(nlam):
        c   = solve_triangular(R[k], QdTMy[k])
        spl = TBSpline.fromTck(xdata,y,(t,np.r_[c,np.zeros(4)],3),edata=edata,natural=True,lam=lams[k],peff=peffs[k])
        spl.logEvidence = logEv[k]
        splinesLam.append(spl)
        ylam.append(Bspl @ c)
        cbar += w[k]*c
    ylam     = np.array(ylam)
    peffMean = w @ peffs
    spline0  = TBSpline.fromTck(xdata,y,(t,np.r_[cbar,np.zeros(4)],3),edata=edata,natural=True,lam='average',peff=peffMean)
    ybar     = w @ ylam
    xmaxs    = xspl[np.argmax(ylam,axis=1)]

    avg = {}
    avg['lams']       = lams
    avg['t']          = t
    avg['M']          = M
    avg['Qd']         = Qd
    avg['R']          = R
    avg['occam']      = occam
    avg['peffs']      = peffs
    avg['peffMean']   = peffMean
    avg['Bspl']       = Bspl
    avg['weights']    = w
    avg['logEv']      = logEv
    avg['splinesLam'] = splinesLam
    avg['splineMean'] = spline0
    avg['ybar']       = ybar
    avg['ylamSpread'] = np.sqrt(w @ (ylam - ybar)**2)
    avg['xmaxSpread'] = np.sqrt(w @ (xmaxs - w @ xmaxs)**2)
    return avg


def _checkAverageArgs(num_knots, fixedKnots, order):
    if (num_knots is not None) or (fixedKnots is not None):
        logger.TBRaise('A smoothing spline has a knot at every unique x. Do not pass num_knots or fixedKnots.')
    if order != 3:
        logger.TBRaise('Smoothing splines are only implemented for order=3.')


def _bootSplineAverage(xdata, ydata, edata, std, useCovariance, Covariance, xspl, num_knots, fixedKnots, 
                       order, numb_samples, seed, nproc, nlam) -> dict:
    """
    bootSpline with lam='average'.
    """
    _checkAverageArgs(num_knots, fixedKnots, order)
    y     = np.asarray(ydata,dtype=float)
    ndata = len(y)
    avg   = _lamAverage(xdata, y, edata, xspl, nlam)
    t     = avg['t']
    ncoeff = len(t)

    fitArgs = {'M':avg['M'],'Qd':avg['Qd'],'R':avg['R'],'occam':avg['occam'],'Bspl':avg['Bspl'],'xspl':xspl}
    samples, mean, err = bootstr_from_gauss(_averageSample,data=[y],data_std_dev=[std],
                                            numb_samples=numb_samples,return_sample=True,seed=seed,
                                            args=fitArgs,nproc=nproc,useCovariance=useCovariance,
                                            Covariance=Covariance)

    # Rebuild the averaged bootstrap splines without refitting.
    splines = [] 
    for sample in samples:
        yBS = sample[:ndata]
        c   = sample[ndata:ndata+ncoeff]
        splines.append(TBSpline.fromTck(xdata,yBS,(t,c,3),edata=edata,natural=True,lam='average',peff=avg['peffMean']))

    res = {}
    res['splineBS']     = splines # in case you want spline functions at bootstrap level
    res['xspl']         = xspl
    res['yspl']         = mean[ndata+ncoeff:-1]
    res['ysple']        = err[ndata+ncoeff:-1]
    res['ysple_lam']    = avg['ylamSpread'] # diagnostic only; already part of ysple
    res['splineMean']   = avg['splineMean']
    # The last entry of each sample is the xmax of the averaged spline.
    res['xmax']         = mean[-1]
    res['xmaxe']        = err[-1]
    res['xmaxe_lam']    = avg['xmaxSpread'] # diagnostic only; already part of xmaxe
    res['lams']         = avg['lams']
    res['peffs']        = avg['peffs']
    res['weights']      = avg['weights']
    res['logEvidences'] = avg['logEv']
    res['splinesLam']   = avg['splinesLam']
    return res


def posteriorSpline(xdata, ydata, edata, numb_samples=300, nsupport=301, seed=DEFAULTSEED, nlam=100) -> dict:
    """
    Smoothing spline with the smoothing parameter lam integrated out, and its Bayesian error band. The
    smoothing spline at fixed lam (see getSpline) is the posterior mean of a Bayesian model, in which 
    the curvature penalty is a Gaussian prior on the B-spline coefficients. Their posterior is Gaussian
    with covariance (B^T W B + lam Omega)^-1. We average over nlam values of lam, uniform in log lam, 
    weighted by their evidence, i.e. integrate lam out with a flat prior in log lam. This is not a 
    bootstrap: the data are used once, as measured, and the error band says how uncertain the true 
    curve is, given the data and the assumption that it is smooth. Unlike the bootstrap error of 
    bootSpline, this includes the uncertainty coming from the smoothing itself.

    Args:
        xdata (array-like)
        ydata (array-like)
        edata (array-like): 
            Errors of ydata, or their covariance matrix.
        numb_samples (int, optional):
            Number of samples drawn from the posterior. Defaults to 300.
        nsupport (int, optional):
            Number of points of the grid xspl, on which the spline is evaluated. Defaults to 301.
        seed (int, optional):
            Seed for the samples. Defaults to DEFAULTSEED.
        nlam (int, optional):
            Number of values of lam. Defaults to 100.

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
    edata = np.asarray(edata,dtype=float)
    y     = np.asarray(ydata,dtype=float)
    xspl  = np.linspace(np.min(xdata),np.max(xdata),nsupport)
    avg   = _lamAverage(xdata, y, edata, xspl, nlam)
    w     = avg['weights']
    t     = avg['t']
    Bspl  = avg['Bspl']
    nB    = avg['R'].shape[2]

    # At fixed lam, the posterior covariance of the coefficients is A^-1 = R^-1 R^-T, so c_k + R^-1 z with 
    # standard normal z is a draw from the posterior. Also accumulate the exact posterior variance on xspl.
    Rinv  = []
    ystat = np.zeros(len(xspl))
    for k in range(nlam):
        Rinvk = solve_triangular(avg['R'][k], np.eye(nB))
        Rinv.append(Rinvk)
        avg['splinesLam'][k].covPost = Rinvk @ Rinvk.T
        ystat += w[k]*np.sum((Bspl @ Rinvk)**2,axis=1)
    ystat = np.sqrt(ystat)

    # Draw lam_k with probability w_k, then the coefficients
    rng = TBRNG(seed)
    ks  = rng.choice(nlam, size=numb_samples, p=w)
    z   = rng.standard_normal((numb_samples, nB))
    splines, xmaxs = [], []
    for i in range(numb_samples):
        k = ks[i]
        c = avg['splinesLam'][k].get_coeffs()[:nB] + Rinv[k] @ z[i]
        splines.append(TBSpline.fromTck(xdata,y,(t,np.r_[c,np.zeros(4)],3),edata=edata,natural=True,lam='average',
                                        peff=avg['peffMean']))
        xmaxs.append(xspl[np.argmax(Bspl @ c)])
    xmaxs = np.array(xmaxs)

    res = {}
    res['xspl']          = xspl
    res['yspl']          = avg['ybar']
    res['ysple_stat']    = ystat
    res['ysple_sys']     = avg['ylamSpread']
    res['ysple']         = np.sqrt(ystat**2 + avg['ylamSpread']**2)
    res['splineMean']    = avg['splineMean']
    res['splineSamples'] = splines
    res['xmax']          = std_median(xmaxs)
    res['xmaxe']         = dev_by_dist(xmaxs)
    res['xmaxe_sys']     = avg['xmaxSpread']
    res['xmaxe_stat']    = np.sqrt(max(res['xmaxe']**2 - avg['xmaxSpread']**2, 0))
    res['lams']          = avg['lams']
    res['peffs']         = avg['peffs']
    res['weights']       = w
    res['logEvidences']  = avg['logEv']
    res['splinesLam']    = avg['splinesLam']
    return res
