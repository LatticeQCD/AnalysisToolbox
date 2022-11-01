# 
# num_int.py                                                               
# 
# D. Clarke
# 
# Some wrappers to do numerical integration.
#

import numpy as np
import scipy.integrate as integrate
import latqcdtools.base.logger as logger
from latqcdtools.math.spline import getSpline


def integrateData(xdata,ydata,method='spline',getQ=False):
    """ Wrapper to integrate data. The default 'spline' will try fitting the data with a cubic spline, then integrate
        the spline function using Gaussian quadrature. """
    try:
        xdata=np.array(xdata)
        ydata=np.array(ydata)
    except:
        logger.TBError("integrateData must be passed array-like objects.")

    if method=='simpson':
        return integrate.simpson(ydata,xdata)

    elif method=='trapezoid':
        return integrate.trapezoid(ydata, xdata)

    elif method=='spline':
        nknots = int(len(xdata)/3)
        data_spline = getSpline(xdata, ydata, num_knots=nknots, rand=False)
        q = data_spline.get_residual()
        if getQ:
            return integrate.quad(data_spline, xdata[0], xdata[-1])[0], q
        else:
            return integrate.quad(data_spline, xdata[0], xdata[-1])[0]

    else:
        logger.TBError("Unknown integration method",method)
