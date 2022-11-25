# 
# math.py                                                               
# 
# D. Clarke
# 
# A handful of wrappers for math functions, as well as functions that I couldn't find in numpy or scipy.
# 
import numpy as np
from scipy.special import kn
import latqcdtools.base.logger as logger
from latqcdtools.base.check import UnderflowError


def fallFactorial(n,m):
    """ Falling factorial n fall to m. """
    if m==0:
        return 1
    if m>n:
        logger.TBError("m>n.")
    prod=1
    for i in range(m):
        prod *= n-i
    return prod


def riseFactorial(n,m):
    """ Rising factorial n rise to m. """
    if m==0:
        return 1
    if n>m:
        logger.TBError("n>m.")
    prod=1
    for i in range(m):
        prod *= n+i
    return prod


# Some methods for dealing with underflow, which is when a number gets so small, it is smaller than the smallest number
# that Python can represent. Functions that drop expoentially are especially prone to underflow. If you like, you can
# protect against underflow by treating these as zero.


def underflowKn(n,x):
    """ Modified Bessel function of the second kind that replaces underflows with 0. """
    try:
        return kn(n,x)
    except UnderflowError:
        return 0.


def underflowExp(x):
    """ Exponential that replaces underflows with 0. """
    try:
        return np.exp(x)
    except UnderflowError:
        return 0.