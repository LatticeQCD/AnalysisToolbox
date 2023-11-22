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
from latqcdtools.base.utilities import isArrayLike 


id_3 = np.eye(3,dtype=complex)
id_4 = np.eye(4,dtype=complex)

ze_3 = np.zeros((3,3), dtype=complex)


def normalize(arr):
    return arr/np.sum(np.abs(arr))


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
