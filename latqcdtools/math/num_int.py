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
from latqcdtools.base.check import checkEqualLengths, checkType
from latqcdtools.base.utilities import envector,unvector


def solveIVP(dydt,t0,tf,y0,method='RK45',args=(),epsrel=1.49e-8,epsabs=1.49e-8) -> np.ndarray:
    """ Wrapper to solve an initial value problem of the form

    dy/dt = dydt(t, y)
    y(t0) = y0

    Args:
        dydt (func): RHS of IVP. Must have signature dydt(t,y).
        t0 (float): initial time 
        tf (float): final time 
        y0 (array-like): y(t0)
        method (str, optional): Integration method. Defaults to 'RK45'.
        args (tuple, optional): Arguments to dydt. Defaults to ().
        epsrel (float, optional): Relative error tolerance. Defaults to 1.49e-8.
        epsabs (float, optional): Absolute error tolerance. Defaults to 1.49e-8.

    Returns:
        array-like: y(tf) 
    """
    a, b = envector(t0, tf)
    checkEqualLengths(a,b)
    def g(A,B):
        sol = integrate.solve_ivp(dydt,(A,B),envector(y0),method=method,args=args,rtol=epsrel,atol=epsabs)
        return unvector(sol.y)[-1]
    h = np.vectorize(g)
    return np.asarray(h(a,b))


def integrateData(xdata,ydata,method='trapezoid'):
    """ Wrapper to integrate data. 

    Args:
        xdata (array-like)
        ydata (array-like)
        method (str, optional): Integration method. Defaults to 'trapezoid'. Possibilities include:
                                > 'simpson' : Simpson rule.
                                > 'trapezoid' : Trapezoidal rule.

    Returns:
        float: Area under ydata. 
    """
    checkType(xdata,'array')
    checkType(ydata,'array')
    xdata=np.array(xdata)
    ydata=np.array(ydata)

    if method=='simpson':
        return integrate.simpson(ydata,xdata)

    elif method=='trapezoid':
        return integrate.trapezoid(ydata, xdata)

    else:
        logger.TBError("Unknown integration method",method)


persistentMethods = ['quad', 'trapezoid', 'romberg']


def integrateFunction(func,a,b,method='persistent',args=(),stepsize=None,limit=1000,epsrel=1.49e-8,epsabs=1.49e-8):
    """ Wrapper to integrate functions. Allows to conveniently adjust the stepsize, and can vectorize scipy.quad, and
        scipy.romberg, which otherwise do not like to handle numpy arrays.

    Args:
        func (func): Integrand. 
        a (float): Lower integration limit.
        b (float): Upper integration limit.
        method (str,optional): Integration method. Defaults to 'persistent_quad_trap'. Possibilities are:
                                > 'persistent' : Try various methods until something works. 
                                > 'quad' : Gaussian quadrature.
                                > 'romberg' : Romberg method. 
                                > 'trapezoid' : Trapezoidal rule.
        args (tuple, optional): Arguments to func. Defaults to ().
        stepsize (float, optional): _description_. Defaults to None.
        limit (int, optional): Upper bound on number of subintervals used in the adaptive algorithm. Defaults to 1000.
        epsrel (float, optional): Relative error tolerance. Defaults to 1.49e-8.
        epsabs (float, optional): Absolute error tolerance. Defaults to 1.49e-8.

    Returns:
        float: Integral of func from a to b. 
    """
    a, b = envector(a,b)
    checkEqualLengths(a,b)
    isVec = False
    if len(a) > 1:
        isVec=True

    if method=='persistent':
        for persistentMethod in persistentMethods:
            try:
                return integrateFunction(func,a,b,args=args,method=persistentMethod,stepsize=stepsize,epsrel=epsrel,epsabs=epsabs)
            except integrate.IntegrationWarning:
                continue
        logger.TBError('No persistent method worked.')

    elif method=='vec_quad':
        def g(A,B):
            return integrate.quad(func, A, B, args=args, limit=limit, epsrel=epsrel, epsabs=epsabs)[0]
        h = np.vectorize(g)
        return np.asarray(h(a,b))

    elif method=='vec_romberg':
        def g(A,B):
            return integrate.romberg(func, A, B, args=args, rtol=epsrel, tol=epsabs)
        h = np.vectorize(g)
        return np.asarray(h(a,b))

    elif method=='quad':
        if isVec:
            return integrateFunction(func,a,b,args=args,method='vec_quad',limit=limit,epsrel=epsrel,epsabs=epsabs)
        else:
            return integrate.quad(func, a, b, args=args, limit=limit, epsrel=epsrel, epsabs=epsabs)[0]

    elif method=='romberg':
        if isVec:
            return integrateFunction(func,a,b,args=args,method='vec_romberg',limit=limit,epsrel=epsrel,epsabs=epsabs)
        else:
            return integrate.romberg(func, a, b, args=args, tol=epsabs, rtol=epsrel) 

    elif method=='trapezoid':
        for i in range(len(b)):
            if b[i] == np.inf or a[i] == -np.inf:
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