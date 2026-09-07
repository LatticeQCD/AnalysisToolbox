#
# pade.py
#
# D. Clarke
#
# Padé approximants. Given the Taylor coefficients of a function about a single point,
# construct a rational approximation that matches the series to as high an order as possible.
#

import numpy as np
from scipy.interpolate import pade
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger
from latqcdtools.math.polynomials import Rational


class singlePointPade:

    """
    Single-point Padé approximant. Given the Taylor coefficients c of some function f about
    a point x0, i.e.
        f(x) = sum_k c[k] (x-x0)**k,
    build the rational function
        R(x) = P(x-x0) / Q(x-x0),   deg(P) = p,  deg(Q) = q,
    whose own Taylor expansion about x0 agrees with f through order p+q. The denominator is
    normalized so that Q(0) = 1. The object is callable and is backed by a Rational.

    The numerator degree p and denominator degree q must both be given. They have to satisfy
    p >= 0, q >= 0, and p + q <= len(c) - 1, since matching the series through order p+q needs
    the coefficients c[0..p+q].

    Example:
        from math import factorial
        c = [1/factorial(k) for k in range(7)]   # Taylor coeffs of exp about 0
        R = singlePointPade(c, p=3, q=3)
        R(0.5)                                    # ~ exp(0.5)
    """

    def __init__(self, c, p, q, x0=0.):
        checkType("array",c=c)
        checkType("int",p=p)
        checkType("int",q=q)
        checkType("real",x0=x0)
        c = np.asarray(c, dtype=float)
        if p < 0 or q < 0:
            logger.TBRaise(f"Need p>=0 and q>=0. Got p={p}, q={q}.")
        if p + q + 1 > len(c):
            logger.TBRaise(f"Taylor series has order {len(c)-1}, which is too low for a [{p}/{q}] "
                           f"Padé approximant. Need len(c) >= p+q+1 = {p+q+1}, but len(c) = {len(c)}.")

        # scipy.interpolate.pade(an, m, n): m is the DENOMINATOR order, n the NUMERATOR order.
        P, Q = pade(c[:p+q+1], q, p)

        # numpy.poly1d stores coefficients highest-power-first; Rational wants ascending order.
        self.num_coeffs = np.asarray(P.coef[::-1], dtype=float)
        self.den_coeffs = np.asarray(Q.coef[::-1], dtype=float)
        self.p        = p
        self.q        = q
        self.x0       = x0
        self.rational = Rational(self.num_coeffs, self.den_coeffs)

    def __repr__(self) -> str:
        # getattr guards: the toolbox logger stringifies self when a TBRaise fires inside __init__,
        # which can happen before these attributes are assigned.
        return "singlePointPade(p={}, q={}, x0={})".format(
            getattr(self,'p','?'), getattr(self,'q','?'), getattr(self,'x0','?'))

    def __call__(self, x):
        scalar = np.isscalar(x)
        y = self.rational(np.asarray(x, dtype=float) - self.x0)
        return float(y) if scalar else y
