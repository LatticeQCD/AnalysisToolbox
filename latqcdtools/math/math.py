# 
# math.py                                                               
# 
# D. Clarke
# 
# A handful of wrappers for math functions, as well as functions that I couldn't find in numpy or scipy.
# 
import numpy as np
from scipy.special import poch
import latqcdtools.base.logger as logger
from latqcdtools.base.check import UnderflowError


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
        return underflowExp(exponent)