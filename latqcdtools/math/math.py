# 
# math.py                                                               
# 
# D. Clarke
# 
# A handful of wrappers for math functions, as well as functions that I couldn't find in numpy or scipy.
# 
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.check import UnderflowError


def fallFactorial(n,m):
    """ Falling factorial n fall to m. """
    if m==0:
        return 1
    if m>n:
        logger.TBError("m>n in falling factorial.")
    prod=1
    for i in range(m):
        prod *= n-i
    return prod


def riseFactorial(n,m):
    """ Rising factorial n rise to m. """
    if m==0:
        return 1
    if n>m:
        logger.TBError("n>m in rising factorial.")
    prod=1
    for i in range(m):
        prod *= n+i
    return prod


def underflowExp(x):
    """ Exponential that replaces underflows with 0. """
    try:
        return np.exp(x)
    except UnderflowError:
        return 0.