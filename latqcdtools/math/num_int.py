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
from latqcdtools.base.utilities import envector
from latqcdtools.math.spline import getSpline


def integrateData(xdata,ydata,method='spline'):
    """ Wrapper to integrate data. The default 'spline' will try fitting the data with a cubic spline, then integrate
        the spline function using Gaussian quadrature. """
    try:
        xdata=np.array(xdata)
        ydata=np.array(ydata)
    except:
        logger.TBError("Must be passed array-like objects.")

    if method=='simpson':
        return integrate.simpson(ydata,xdata)

    elif method=='trapezoid':
        return integrate.trapezoid(ydata, xdata)

    elif method=='spline':
        nknots = int(len(xdata)/3)
        data_spline = getSpline(xdata, ydata, num_knots=nknots, rand=False)
        return integrate.quad(data_spline, xdata[0], xdata[-1])[0]

    else:
        logger.TBError("Unknown integration method",method)


def integrateFunction(func,a,b,method='persistent_quad_trap',args=(),stepsize=None,limit=1000,epsrel=1e-8):
    """ Wrapper to integrate functions. Allows to conveniently adjust the stepsize, and can vectorize scipy.quad,
        which otherwise does not like to handle numpy arrays. """
    a = envector(a)
    b = envector(b)
    if not len(a) == len(b):
        logger.TBError('Integration bounds must have the same length.')

    if method=='persistent_quad_trap':
        try:
            return integrateFunction(func,a,b,args=args,method='quad',stepsize=None)
        except integrate.IntegrationWarning:
            return integrateFunction(func,a,b,args=args,method='trapezoid',stepsize=None)

    elif method=='vec_quad':
        def g(A,B):
            # TODO: put a check to make sure the integration error is small compared to result
            return integrate.quad(func, A, B, args=args, limit=limit, epsrel=epsrel)[0]
        h = np.vectorize(g)
        if len(a)==1:
            return np.asarray(h(a,b))[0]
        else:
            return np.asarray(h(a,b))

    elif method=='quad':
        return integrate.quad(func, a, b, args=args, limit=limit, epsrel=epsrel)[0]

    elif method=='trapezoid':
        for i in range(len(b)):
            if b[i] == np.inf:
                logger.TBError('Trapezoid rule is meant for definite integrals.')
        if stepsize is None:
            x = np.array([ np.linspace(a[i], b[i], 101) for i in range(len(b)) ])
        else:
            x = np.array([ np.arange(start=a[i], stop=b[i], step=stepsize) for i in range(len(b)) ])
        y = func(x,*args)
        if len(a)==1:
            return integrateData(x, y, method='trapezoid')[0]
        else:
            return integrateData(x, y, method='trapezoid')

    else:
        logger.TBError('Unrecognized integration method',method)