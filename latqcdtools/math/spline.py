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
    """ Return a list of nknots evenly spaced knots. """
    if len(xdata)<nknots:
        logger.TBError('number of data < number of knots')
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
    """ Return a list of nknots randomly spaced knots. """
    rng = np.random.default_rng(SEED)
    flat_xdata = np.sort(np.asarray(xdata))
    sample_xdata = rng.choice(flat_xdata,int(nknots+1+(1-randomization_factor)*(len(flat_xdata)-nknots)),
                              replace=False)
    # Retry if too many data points are removed by np.unique
    if len(np.unique(sample_xdata)) < nknots + 1:
        return _random_knots(xdata, nknots, randomization_factor)
    return _even_knots(sample_xdata, nknots)


class customSpline:

    def __init__(self,xdata,ydata,edata=None,knots=None,order=3,smooth=None):
        if edata is None:
            weights = None
        else:
            weights = 1/edata
            weights = np.r_[weights[0],weights,weights[-1]]
        # This is a trick to try to enforce no curvature at the endpoints.
        xx = np.r_[xdata[0]  ,xdata  ,xdata[-1]]
        yy = np.r_[ydata[0]  ,ydata  ,ydata[-1]] 
        self.tck = splrep(xx, yy, t=knots, k=order, w=weights, s=smooth)

    def __call__(self,x):
        return splev(x,self.tck)

    def get_coeffs(self):
        return self.tck[1]


def getSpline(xdata, ydata, num_knots=None, edata=None, order=3, rand=False, fixedKnots=None, 
              getAICc=False, natural=False, smooth=None):
    """ This is a wrapper that calls SciPy spline fitting methods, depending on your needs. Calls LSQUnivariateSpline
    by default. If you need to ensure a well defined second derivative at the knots, we call instead UnivariateSpline,
    since LSQUnivariate spline seems to have no smoothing option. Sadly if you call UnivariateSpline, you can't specify
    the knots, so you can't both pick knots and smooth.

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
        logger.TBError('len(xdata), len(ydata) =',len(xdata),len(ydata))

    if natural:
        if num_knots is not None:
            logger.TBError('Scipy natural spline chooses knots automatically.')
        if order != 3:
            logger.TBError("Scipy natural spline only implemented for cubic splines.")
        spline = CubicSpline(x=xdata,y=ydata,bc_type='natural')

    else:
        checkType(num_knots,int)
        nknots = num_knots
        if fixedKnots is not None:
            if type(fixedKnots) is not list:
                logger.TBError("knots must be specified as a list.")
            nknots -= len(fixedKnots)
            if nknots < 0:
                logger.TBError("len(fixedKnots)",len(fixedKnots),"exceeds num_knots",num_knots)
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
            logger.TBError("You can't put a knot to the left of the x-data. knots, xdata[0] = ",knots,xdata[0])
        if knots[-1]>xdata[-1]:
            logger.TBError("You can't put a knot to the right of the x-data. knots, xdata[-1] = ",knots,xdata[-1])
        spline = customSpline(xdata, ydata, edata=edata, knots=knots, order=order, smooth=smooth)

    if getAICc:
        cov = np.diag(edata**2)
        return spline, AICc(xdata, ydata, cov, spline)
    else:
        return spline


def getSplineErr(xdata, xspline, ydata, ydatae, num_knots=None, order=3, rand=False, fixedKnots=None, natural=False):
    """ Use getSpline to smooth mean and error bars. Create a spline-smooth band from that. """
    spline_lower = getSpline(xdata, ydata - ydatae, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural)(xspline)
    spline_upper = getSpline(xdata, ydata + ydatae, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural)(xspline)
    spline_center = getSpline(xdata, ydata, num_knots=num_knots, order=order, rand=rand, fixedKnots=fixedKnots, natural=natural)(xspline)
    spline_err = ((spline_upper - spline_center) + (spline_center - spline_lower)) / 2
    return spline_center, spline_err

