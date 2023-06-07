# 
# math.py                                                               
# 
# D. Clarke
# 
# Wrappers for some math functions, some math functions that I couldn't find in numpy or scipy, and some methods
# for comparing math objects. 
#


import numpy as np
from scipy.special import poch
import latqcdtools.base.logger as logger
from latqcdtools.base.check import UnderflowError
from latqcdtools.base.utilities import envector, isArrayLike 


def fallFactorial(n,m):
    """ Falling factorial n fall to m. """
    if m>n:
        logger.TBError("m>n.")
    return poch(n-m+1,m)


def riseFactorial(n,m):
    """ Rising factorial n rise to m. """
    if n>m:
        logger.TBError("n>m.")
    return poch(n,m)


def logDet(mat):
    """ Logarithm of determinant. """
    _, ans = np.linalg.slogdet(mat)
    return ans


# Some methods for dealing with underflow, which is when a number gets so small, it is smaller than the smallest number
# that Python can represent. If you like, you can suppress underflow errors with these functions.


def underflowExp(x):
    """ Exponential that replaces underflows with 0. """
    try:
        return np.exp(x)
    except UnderflowError:
        # log(2**(-1022)) ~ -708
        flowMask = x>-708
        return flowMask*np.exp(flowMask*x)


def underflowPower(x,n):
    """ x**n that replaces underflows with 0. """
    try:
        return np.power(x,n)
    except UnderflowError:
        return np.power(np.log(np.exp(x)),n)


def underflowMultiply(*args):
    """ Product of arbitrary number of reals/np.arrays that replaces underflows with 0. """
    prod = 1
    try:
        for x in args:
            prod *= x
        return prod
    except UnderflowError:
        exponent = 0
        for x in args:
            exponent += np.log(x)
        return underflowExp(exponent)#


# Some methods for doing comparisons between math objects.


def rel_check(a, b, prec = 1e-6, abs_prec = 1e-14):
    """ Check whether a and b are equal. a and b can be array-like, float-like, or complexes. If a
    and b are array-like, we check that they are element-wise equal within the tolerance. 

    Args:
        a (obj)
        b (obj)
        prec (float, optional): Relative precision. Defaults to 1e-6.
        abs_prec (float, optional): Absolute precision. Defaults to 1e-14.

    Returns:
        bool: True if a and b are equal. 
    """
    if isArrayLike(a):
        if np.shape(a) != np.shape(b):
            logger.TBError('a and b must have the same shape. Received a, b shapes =',np.shape(a),np.shape(b))
        return np.allclose( a, b, rtol = prec, atol = abs_prec)
    else:
        try:
            return np.isclose( a, b, rtol = prec, atol = abs_prec)
        except TypeError:
            logger.TBError('Expected reals, complexes, or array-like. Received a, b types =',type(a),',',type(b))


def print_results(res, res_true, res_err = None, res_err_true = None, text = "", prec = 1e-10, abs_prec = None):
    """ Compares element-by-element the results of res with res_true. (Does the same with res_err and res_err_true,
        if you like.) Carries out with precision prec. Use abs_prec for comparisons with zero. """
    test = True
    compareZero = False

    if abs_prec is not None:
        compareZero = True
    else:
        abs_prec = 1e-14

    res = envector(res)
    res_true = envector(res_true)
    if res_err is not None:
        res_err = envector(res_err)
    if res_err_true is not None:
        res_err_true = envector(res_err_true)

    for i in range(len(res)):
        if not rel_check(res[i], res_true[i], prec, abs_prec):
            test = False
            logger.info("res[" + str(i) + "] = " + str(res[i]) + " != res_true[" + str(i) + "] = " + str(res_true[i]))
        if res_err is not None and res_err_true is not None:
            if not rel_check(res_err[i], res_err_true[i], prec, abs_prec):
                test = False
                logger.info("res_err[" + str(i) + "] = " + str(res_err[i]) + " != res_err_true[" + str(i) + "] = " + str(res_err_true[i]))

    if test:
        if compareZero:
            logger.TBPass(text,'(abs_prec = %.2e)' % abs_prec)
        else:
            logger.TBPass(text,'(prec = %.2e)' % prec)
    else:
        if compareZero:
            logger.TBFail(text, '(abs_prec = %.2e)' % abs_prec)
        else:
            logger.TBFail(text, '(prec = %.2e)' % prec)
