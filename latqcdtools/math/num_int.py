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


# TODO: implement romberg?
def integrateFunction(func,a,b,method='persistent_quad_trap',stepsize=None):
    """ Wrapper to integrate functions. """
    if method=='persistent_quad_trap':
        try:
            return integrate.quad(func, a, b, limit=1000, epsrel=1e-8)[0]
        except integrate.IntegrationWarning:
            return integrateFunction(func,a,b,method='trapezoid',stepsize=None)
        # TODO: can you vectorize quad?
    elif method=='quad':
        return integrate.quad(func, a, b, limit=1000, epsrel=1e-8)[0]
    elif method=='trapezoid':
        # TODO: try changing b until you hit some relative tolerance? Take OOM steps?
        if b == np.inf:
            b = 10 * a
        if stepsize is None:
            x = np.linspace(a, b, 101)
        else:
            x = np.arange(start=a, stop=b, step=stepsize)
        y = func(x)
        return integrateData(x, y, method='trapezoid')

#        import numpy as np
#        import scipy.integrate
#
#        def new_quad(f, a, b):
#            def g(a, b):
#                return scipy.integrate.quad(f, a, b)
#
#            h = np.vectorize(g)
#            res = h(a, b)
#            # np.asarray ensures that we don't return 0d arrays
#            # tuple is in case we care about returning a tuple
#            return tuple(np.asarray(res))
#
#        def f(x):
#            return x
#
#        a = 0
#        b = [1, 2, 3]
#        new_quad(f, a, b)  # (array([0.5, 2. , 4.5]), array([5.55111512e-15, 2.22044605e-14, 4.99600361e-14]))