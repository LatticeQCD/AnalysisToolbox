# 
# spline.py                                                               
# 
# D. Clarke, H. Sandmeyer
# 
# Generally speaking, one should use scipy's methods for splines, like interp1d, UnivariateSpline, etc. However
# it is a bit inconvenient to use when one wants control over the knots and endpoints. That is what this module is for.
# 
import numpy as np
from scipy.interpolate import LSQUnivariateSpline
import latqcdtools.base.logger as logger

def even_knots(xdata, nknots):
    """ Return a list of nknots evenly spaced knots. """
    try:
        flat_xdata = np.sort(np.concatenate(xdata))
    except ValueError:
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


def random_knots(xdata, nknots, randomization_factor=1, SEED=None):
    """ Return a list of nknots randomly spaced knots. """
    np.random.seed(SEED)
    try:
        flat_xdata = np.sort(np.concatenate(xdata))
    except ValueError:
        flat_xdata = np.asarray(xdata)
    sample_xdata = np.random.choice(flat_xdata,int(nknots+1+(1-randomization_factor)*(len(flat_xdata)-nknots)),
                                    replace=False)
    # Retry if too many data points are removed by np.unique
    if len(np.unique(sample_xdata)) < nknots + 1:
        return random_knots(xdata, nknots, randomization_factor)
    ret = even_knots(sample_xdata, nknots)
    return ret


def getSpline(xdata, ydata, num_knots, order=3, rand=False, fixedKnots=None):
    if type(num_knots) is not int:
        logger.TBError("Please specify an integer number of knots.")
    nknots = num_knots
    if fixedKnots is not None:
        if type(fixedKnots) is not list:
            logger.TBError("knots must be specified as a list.")
        nknots -= len(fixedKnots)
        if nknots < 0:
            logger.TBError("len(fixedKnots) cannot exceed num_knots.")
    if nknots>0:
        if rand:
            knots = random_knots(xdata,nknots)
        else:
            knots = even_knots(xdata,nknots)
    else:
        knots=[]
    if fixedKnots is not None:
        for knot in fixedKnots:
            knots.append(knot)
    knots = sorted(knots)
    if knots[0]<xdata[0]:
        logger.TBError("You can't but a knot to the left of the x-data...")
    logger.debug("Knots are:",knots)
    spline = LSQUnivariateSpline(xdata, ydata, knots, k=order)
    return spline