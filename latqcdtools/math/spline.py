# 
# spline.py                                                               
# 
# D. Clarke
# 
# Generally speaking, one should use scipy's methods for splines, like interp1d, UnivariateSpline, etc. However
# it is a bit inconvenient to use when one wants control over the knots and endpoints. That is what this module is for.
# 

import numpy as np
from scipy.interpolate import CubicSpline, splrep, splev
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import AICc
from latqcdtools.base.check import checkType


def _even_knots(xdata, nknots):
    """ 
    Return a list of nknots evenly spaced knots. 
    """
    if len(xdata)<nknots:
        logger.TBRaise('number of data < number of knots')
    flat_xdata = np.sort(np.asarray(xdata))
    # to ensure no knot sits at a data position
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
    Return a list of nknots randomly spaced knots. 
    """
    rng = np.random.default_rng(SEED)
    flat_xdata = np.sort(np.asarray(xdata))
    sample_xdata = rng.choice(flat_xdata,int(nknots+1+(1-randomization_factor)*(len(flat_xdata)-nknots)),
                              replace=False)
    # Retry if too many data points are removed by np.unique
    if len(np.unique(sample_xdata)) < nknots + 1:
        return _random_knots(xdata, nknots, randomization_factor)
    return _even_knots(sample_xdata, nknots)


class TBSpline:

    """
    A class that prepares a splrep and wraps it with splev.
    """

    def __init__(self,xdata,ydata,edata=None,knots=None,order=3,naturalLike=False):

        self.xspline = np.copy(xdata)
        self.yspline = np.copy(ydata)

        # A "natural spline" is technically a solve that enforces zero curvature at the
        # endpoints. This class on the other hand wraps smoothing splines, usually being
        # applied to data with errors. A strategy to impose zero curvature at the endpoints
        # is to create a fake datum before each endpoint such that the datum is colinear
        # with both the endpoint and the next-innermost point. When there are weights,
        # we weight these last three points more than the rest of the data to try to force
        # the spline to pass through them.
        if naturalLike:
            dxL = xdata[1 ]-xdata[0 ]
            dxR = xdata[-1]-xdata[-2]
            dyL = ydata[1 ]-ydata[0 ]
            dyR = ydata[-1]-ydata[-2]
            self.xspline = np.r_[xdata[0]-dxL, xdata, xdata[-1]+dxR]
            self.yspline = np.r_[ydata[0]-dyL, ydata, ydata[-1]+dyR]

        if edata is None:
            self.weights = None
            smooth = None
        else:
            self.weights = 1/edata
            smooth = len(self.weights)
            if naturalLike:
                weightL = self.weights[0]
                weightR = self.weights[-1] 
                self.weights = np.r_[weightL, self.weights, weightR] 

        self.tck = splrep(self.xspline, self.yspline, t=knots, k=order, w=self.weights, s=smooth, task=-1)

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


def getSpline(xdata, ydata, num_knots=None, edata=None, order=3, rand=False, fixedKnots=None, 
              getAICc=False, natural=False):
    """ 
    This is a wrapper that calls SciPy spline interpolation methods, depending on your needs. Generally
    this uses scipy.interpolate.splrep, which uses B-splines. If natural=True and edata=None, it will
    use scipy.interpolate.CubicSpline to solve. If natural=True and edata are provided, it will do a
    smoothing spline that attempts to force no curvature at the endpoints, based on the penultimate points. 

    Args:
        xdata (array-like)
        ydata (array-like)
        num_knots (int):
            The number of knots.
        edata (array-like, optional): 
            Error data. Defaults to None.
        order (int, optional):
            Order of the spline. Defaults to 3.
        rand (bool, optional): 
            Use randomly placed knots? Defaults to False.
        fixedKnots (array-like, optional):
            List of user-specified knots. Defaults to None.
        getAICc (bool, optional): 
            Return corrected Aikake information criterion? Defaults to False.
        natural (bool, optional): 
            Try a natural (no change in slope at the endpoints) cubic spline. Defaults to False. 

    Returns:
        callable spline object
        AICc (optionally)
    """

    if len(xdata) != len(ydata):
        logger.TBRaise('len(xdata), len(ydata) =',len(xdata),len(ydata))
    if natural and (order != 3):
        logger.TBRaise("Natural splines have order=3 by definition.")

    if natural and (edata is None): 
        if num_knots is not None:
            logger.TBRaise('Natural spline without error is a solve that chooses knots automatically.')
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
                knots = _random_knots(xdata,nknots)
            else:
                knots = _even_knots(xdata,nknots)
        else:
            knots=[]
        if fixedKnots is not None:
            for knot in fixedKnots:
                knots.append(knot)
        knots = sorted(knots)
        if knots[0]<xdata[0]:
            logger.TBRaise("You can't put a knot to the left of the x-data. knots, xdata[0] = ",knots,xdata[0])
        if knots[-1]>xdata[-1]:
            logger.TBRaise("You can't put a knot to the right of the x-data. knots, xdata[-1] = ",knots,xdata[-1])
        spline = TBSpline(xdata, ydata, edata=edata, knots=knots, order=order, naturalLike=natural)

    if getAICc:
        cov = np.diag(1/spline.weights**2)
        return spline, AICc(spline.xspline, spline.yspline, cov, spline)
    else:
        return spline


def getSplineErr(xdata, xspline, ydata, ydatae, num_knots=None, order=3, rand=False, fixedKnots=None, natural=False):
    """ 
    Use getSpline to smooth mean and error bars. Create a spline-smooth band from that. 
    """
    spline_lower = getSpline(xdata, ydata - ydatae, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural)(xspline)
    spline_upper = getSpline(xdata, ydata + ydatae, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural)(xspline)
    spline_center = getSpline(xdata, ydata, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural)(xspline)
    spline_err = ((spline_upper - spline_center) + (spline_center - spline_lower)) / 2
    return spline_center, spline_err

