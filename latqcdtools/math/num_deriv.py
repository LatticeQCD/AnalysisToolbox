#
# num_deriv.py
#
# H. Sandmeyer
#
# Some methods for calculating numerical derivatives. All use the central difference strategy, which has O(h^2)
# corrections, where h is the step size.
#

import numpy as np


def _best_h(x):
    """ This routine attempts to automatically determine the best possible step size h when calculating numerical
    derivatives. One does not like the step size to be too large, since one pays a ~h^2 penalty. On the other hand,
    if the step size is too small, the function may not change within machine precision eps, and hence the derivative
    will be incorrectly estimated.

    Let M_0, M_3 > 0 s.t. |f| < M_0 and |f'''| < M_3. Let f_computed = f + f*eps, where eps is the machine precision.
    For central differences,
        ( f(x+h) - f(x-h) )/2h - f' = h^3 f'''(x)/3.
    Using the triangle inequality, one can construct an upper bound of the LHS (when calculated on a computer). The h
    which minimizes this upper bound is
        h = (3 eps M_0/M_3)^(1/3).
    Sadly this requires knowledge of M_3, which we presumably don't have since we resorted to numerical derivatives in
    the first place. A compromise that has some memory of this optimization problem is
        h ~ x eps^(1/3).

    We choose eps = 1.1e-16, which is roughly the difference between 1.0 and the next-smallest representable float less
    than 1.0 in 64-bit. The difference between 1.0 and the next-largest float is slightly bigger. """
    eps   = 1.1e-16
    small = pow(eps,1/3)
    h     = small*(abs(x) + small) # Want something meaningful also when x = 0.
    return h


def diff_deriv(x, func, args = (), h = None):
    """ Numerical derivative using central difference. """
    if h is None:
        h = _best_h(x)
    up   = x + h
    down = x - h
    return (func(up, *args) - func(down, *args)) / (2*h)


def diff_grad(params, func, args = (), h = None) -> np.ndarray:
    """ Gradient using difference quotient. """
    ret = [0.0]*len(params)
    up = np.array(params, dtype = float)
    down = np.array(params, dtype = float)
    for i in range(len(params)):
        if h is None:
            h = _best_h(params[i])
        up[i] += h
        down[i] -= h
        ret[i] = (func(up, *args) - func(down, *args)) / (2*h)
        up[i] = params[i]
        down[i] = params[i]
    return np.array(ret)


def diff_hess(params, func, args = (), h = None) -> np.ndarray:
    """ Hessian using difference quotient. """

    # This has to be a list, as we might put in arrays, if params is higher dimensional
    ret = [ [0.0]*len(params) for _ in range(len(params)) ]
    # convert to float
    up = np.array(params, dtype = float)
    down = np.array(params, dtype = float)
    upup = np.array(params, dtype = float)
    updown = np.array(params, dtype = float)
    downdown = np.array(params, dtype = float)
    downup = np.array(params, dtype = float)

    for i in range(len(params)):
        if h is None:
            hi = _best_h(params[i])
        else:
            hi = h
        up[i] += hi
        upup[i] = up[i]
        updown[i] = up[i]
        down[i] -= hi
        downdown[i] = down[i]
        downup[i] = down[i]
        ret[i][i]=(func(up,*args)+func(down,*args)-2*func(params,*args))/(4*(hi/2.0)**2)

        for j in range(i):
            if h is None:
                hj = _best_h(params[j])
            else:
                hj=h

            upup[j]+=hj
            updown[j]-=hj
            downdown[j]-=hj
            downup[j]+=hj

            ret[i][j]=(func(upup,*args) + func(downdown,*args)-func(updown,*args) - func(downup,*args)) / (4*hi*hj)
            ret[j][i]=ret[i][j]

            upup[j]=params[j]
            updown[j]=params[j]
            downdown[j]=params[j]
            downup[j]=params[j]

        up[i] = params[i]
        upup[i] = up[i]
        updown[i] = up[i]
        down[i] = params[i]
        downdown[i] = down[i]
        downup[i] = down[i]

    return np.array(ret)


def diff_fit_grad(x, params, func, args = (), h = None):
    """ When fitting we're trying to optimize params, and hence we want to think of func as a function of
    its parameters rather than x. """ 
    def f(p):
        return func(x, p, *args)
    return diff_grad(params, f, h = h)


def diff_fit_hess(x, params, func, args = (), h = None):
    """ When fitting we're trying to optimize params, and hence we want to think of func as a function of
    its parameters rather than x. """ 
    def f(p):
        return func(x, p, *args)
    return diff_hess(params, f, h = h)


def diff_jac(params, func, args = (), h = None) -> np.ndarray:
    return diff_grad(params, func, args, h).transpose()
