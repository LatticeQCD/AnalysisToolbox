# 
# spline.py                                                               
# 
# D. Clarke, Claude Code
# 
# Wrappers for SciPy splines that give control over knot placement and endpoints. Depending on the arguments,
# getSpline gives either an interpolating natural cubic spline or a least-squares (regression) spline with
# fixed knots. 
# 

import numpy as np
from scipy.interpolate import CubicSpline, BSpline, splrep, splev
from scipy.linalg import null_space
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import AICc
from latqcdtools.base.check import checkType, checkEqualLengths
from latqcdtools.base.utilities import toNumpy, find_nearest_idx
from latqcdtools.base.initialize import TBRNG, DEFAULTSEED
from latqcdtools.statistics.statistics import std_median, dev_by_dist, countParams
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
    rng = np.random.default_rng(SEED)
    flat_xdata = np.sort(np.asarray(xdata))
    nsample = int(nknots+1+(1-randomization_factor)*(len(flat_xdata)-nknots))
    sample_xdata = rng.choice(flat_xdata,nsample,replace=False)
    # Retry if too many data points are removed by np.unique
    while len(np.unique(sample_xdata)) < nknots + 1:
        sample_xdata = rng.choice(flat_xdata,nsample,replace=False)
    return _even_knots(sample_xdata, nknots)


class TBSpline:

    """
    Least-squares (regression) spline with fixed interior knots. Without natural, prepares a splrep
    with task=-1 and wraps it with splev. With natural, does the least-squares fit itself in the
    subspace of cubic B-splines with zero curvature at the endpoints.
    """

    def __init__(self,xdata,ydata,edata=None,knots=None,order=3,natural=False):

        self.xspline = np.copy(xdata)
        self.yspline = np.copy(ydata)
        self.natural = natural

        if edata is None:
            self.weights = None
            smooth = None
        else:
            self.weights = 1/np.asarray(edata)
            # Ignored by splrep when task=-1.
            smooth = len(self.weights)

        if natural:
            self.tck = self._naturalFit(knots,order)
        else:
            self.tck = splrep(self.xspline, self.yspline, t=knots, k=order, w=self.weights, s=smooth, task=-1)

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

    def get_nparams(self) -> int:
        """
        Number of free parameters. Natural boundary conditions remove 2.
        """
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
              getAICc=False, natural=False, seed=None):
    """ 
    This is a wrapper that calls SciPy spline methods, depending on your needs. If natural=True and
    edata=None, it returns a natural cubic interpolating spline from scipy.interpolate.CubicSpline,
    with a knot at every data point. Otherwise it returns a TBSpline, i.e. a (weighted, if edata are
    given) least-squares B-spline fit with fixed knots. If natural=True and edata are provided, this
    least-squares spline has zero curvature at the endpoints. 

    Args:
        xdata (array-like)
        ydata (array-like)
        num_knots (int):
            The number of interior knots, including fixedKnots. Not used for the interpolating spline.
        edata (array-like, optional): 
            Error data. Defaults to None.
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

    if natural and (edata is None): 
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
    # Both splines must use the same random knots.
    if rand and (seed is None):
        seed = int(TBRNG().integers(2**32))
    spline_lower  = getSpline(xdata, ydata - ydatae, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural, seed=seed)(xspline)
    spline_upper  = getSpline(xdata, ydata + ydatae, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural, seed=seed)(xspline)
    spline_center = (spline_lower+spline_upper)/2 
    spline_err    = (spline_upper-spline_lower)/2 
    return spline_center, spline_err


def bootSpline(xdata, ydata, edata, num_knots=None, order=3, rand=False, fixedKnots=None, 
               natural=False, numb_samples=300, nsupport=301, seed=DEFAULTSEED) -> dict:
    """
    Given xdata, ydata, edata, create a spline. Use a Gaussian (parametric) bootstrap, drawing ydata
    from normal distributions of width edata, to propagate uncertainties of the data into an error band
    for the spline. Gives back a dictionary whose xspl, yspl, and ysple entries can be used to plot a
    spline with error bars. yspl is the median over bootstrap samples and ysple the corresponding
    spread (dev_by_dist). splineMean is the spline fit to the original data, not an average. If
    rand=True, each bootstrap sample gets its own random knots.
    """
    checkType("int",seed=seed)
    checkEqualLengths(xdata,ydata,edata)
    rng     = TBRNG(seed)
    xspl    = np.linspace(np.min(xdata),np.max(xdata),nsupport)
    spline0 = getSpline(xdata=xdata,ydata=ydata,edata=edata,num_knots=num_knots,order=order,
                        rand=rand,fixedKnots=fixedKnots,natural=natural,seed=seed)
    nparams = countParams(spline0,params=())
    ndata   = len(ydata)
    splines = [] 
    AICcs   = []
    xmaxs   = []
    res     = {}

    iboot = 0
    while iboot<numb_samples:
        yBS = rng.normal(ydata,edata)
        if rand:
            knotSeed = int(rng.integers(2**32))
        else:
            knotSeed = None
        if nparams>=ndata:
            spl = getSpline(xdata=xdata,ydata=yBS,num_knots=num_knots,edata=edata,order=order,rand=rand,
                            fixedKnots=fixedKnots,natural=natural,seed=knotSeed)
        else:
            spl, AICc = getSpline(xdata=xdata,ydata=yBS,num_knots=num_knots,edata=edata,order=order,rand=rand,
                                  fixedKnots=fixedKnots,natural=natural,getAICc=True,seed=knotSeed)
            AICcs.append(AICc)
        splines.append(spl)

        max_idx = find_nearest_idx(spl(xspl),np.max(spl(xspl)))
        xmaxs.append(xspl[max_idx])

        iboot += 1

    ys, yes = [], []
    for x in xspl:
        splx = []
        for iboot in range(numb_samples):
            splx.append(splines[iboot](x))
        splx = np.array(splx)
        ys.append(std_median(splx))
        yes.append(dev_by_dist(splx))
    ys, yes, AICcs, xmaxs = toNumpy(ys, yes, AICcs, xmaxs)

    res['AICcs']      = AICcs   # in case you want to diagnose fit quality; empty if nparams >= ndata
    res['splineBS']   = splines # in case you want spline functions at bootstrap level
    res['xspl']       = xspl
    res['yspl']       = ys
    res['ysple']      = yes
    res['splineMean'] = spline0
    res['xmax']       = std_median(xmaxs)  # where is the maximum on xspl grid
    res['xmaxe']      = dev_by_dist(xmaxs)
                        
    return res
