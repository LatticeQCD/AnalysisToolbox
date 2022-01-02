#
# num_deriv.py
#
# H. Sandmeyer
#
# Some methods for calculating numerical derivatives.
#
try: # TODO: remove this and just make it part of the requirements.txt right?
    from algopy import UTPM
    import algopy
    algopy_avail = True
except ImportError:
    algopy_avail = False
import numpy as np


# alg_grad with algopy
def alg_grad(x, func, args = ()):
    if not algopy_avail:
        raise ImportError("Algopy not installed. Please install algopy or use diff_grad")
    x = UTPM.init_jacobian(x)
    y = func(x, *args)
    return UTPM.extract_jacobian(y)


def func_part(x, func, args):
    func_res = func(x, *args)
    z = algopy.zeros(len(func_res), dtype=x)
    for i in range(len(func_res)):
        z[i] = func_res[i]
    return z

# alg_jac using algopy
def alg_jac(x, func, args = ()):
    try:
        func(x, *args)[0]
    except IndexError:
        # If func returns a scalar, the jacobian is just the gradient
        return alg_grad(x, func, args)
    if not algopy_avail:
        raise ImportError("Algopy not installed. Please install algopy or use diff_jac")
    x = UTPM.init_jacobian(x)
    y = func_part(x, func, args)
    return UTPM.extract_jacobian(y)


def best_eps(x):
    num_prec = pow(1.1e-16, 1 / 3.0)
    eps = num_prec * (abs(x) + num_prec)
    return eps


def diff_grad(params, func, args = (), eps = None, expand = False):
    """ Gradient using difference quotient. """
    ret = [ 0.0 for i in range(len(params))]
    up = np.array(params, dtype = float)
    down = np.array(params, dtype = float)
    for i in range(len(params)):
        if eps is None:
            eps = best_eps(params[i])
        up[i] += eps
        down[i] -= eps
        if expand:
            ret[i] = (func(*(tuple(up) + tuple(args)))
                    - func(*(tuple(down) + tuple(args)))) / (2*eps)
        else:
            ret[i] = (func(up, *args) - func(down, *args)) / (2*eps)
        up[i] = params[i]
        down[i] = params[i]
    return np.array(ret)


# Gradient is already Jacobian
def diff_jac(params, func, args = (), eps = None, expand = False):
    return diff_grad(params, func, args, eps, expand).transpose()