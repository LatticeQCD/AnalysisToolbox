# 
# statisticalPhysics.py                                                               
# 
# D. Clarke
# 
# A collection of methods relevant for statistical physics calculations.
#


import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.statistics.statistics import std_mean


def _printExponent(prefix, exponent):
    if exponent is not None:
        logger.info(prefix, round(exponent, 4))


class UniversalityClass:

    """ 
    Skeleton universality class from which all others inherit.
    """

    symm  = None
    d     = None
    alpha = None
    beta  = None
    gamma = None
    delta = None
    nu    = None
    eta   = None
    omega = None

    def __repr__(self) -> str:
        return "UniversalityClass"

    def name(self) -> str:
        return str(self.d)+"d, "+str(self.symm)

    def exponentSummary(self):
        logger.info()
        logger.info("Summary of "+self.name()+" critical exponents:")
        _printExponent(" alpha =",self.alpha)
        _printExponent("  beta =",self.beta)
        _printExponent(" gamma =",self.gamma)
        _printExponent(" delta =",self.delta)
        _printExponent(" omega =",self.omega)
        _printExponent("    nu =",self.nu)
        logger.info()

    def hyperscalingCheck(self, tol=1e-6) -> bool:
        err1 = 2*self.beta+self.gamma-2+self.alpha
        err2 = 2*self.beta*self.delta-self.gamma-2+self.alpha
        err3 = self.nu*self.d-2+self.alpha
        if err1 > tol:
            logger.TBFail(self.name(),"fails hyperscaling check 1, err =",err1)
            return False
        if err2 > tol:
            logger.TBFail(self.name(),"fails hyperscaling check 2. err =",err2)
            return False
        if err3 > tol:
            logger.TBFail(self.name(),"fails hyperscaling check 3. err =",err3)
            return False
        return True


# 3d XY model
class O2_3d(UniversalityClass):
    """ 
    3d O(2) critical exponents from JHEP08 (2016) 036 
    """
    symm  = "O(2)"
    d     = 3
    eta   = 0.03852
    nu    = 0.6719
    gamma = nu*(2-eta)
    beta  = 0.5*(nu*d-gamma)
    alpha = 2 - nu*d
    delta = (gamma+2*alpha)/(2*beta)
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name


# 3d Heisenberg model
class O3_3d(UniversalityClass):
    """ 
    3d O(3) critical exponents from JHEP08 (2016) 036 
    """
    symm  = "O(3)"
    d     = 3
    eta   = 0.0386 
    nu    = 0.7121 
    gamma = nu*(2-eta)
    beta  = 0.5*(nu*d-gamma)
    alpha = 2 - nu*d
    delta = (gamma+2*alpha)/(2*beta)
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name


class O4_3d(UniversalityClass):
    """ 
    3d O(4) critical exponents from Nucl. Phys. B 675, 533-554 (2003). 
    """
    symm  = "O(4)"
    d     = 3
    beta  = 0.380
    delta = 4.824
    alpha = 2.-beta*(1.+delta)
    gamma = beta*(delta-1.)
    nu    = (beta/d)*(1+delta)
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name


# 3d Ising model
class Z2_3d(UniversalityClass):
    """ 
    3d Z_2 critical exponents from J. Stat. Phys. 157. 869-914 (2014). 
    """
    symm  = "Z_2"
    d     = 3
    nu    = 0.62999
    eta   = 0.03631
    alpha = 2 - nu*d
    gamma = nu*(2 - eta)
    beta  = (nu*d -gamma)/2
    delta = nu*d/beta-1.
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name


# 2d Ising model
class Z2_2d(UniversalityClass):
    """ 
    Exact solution for 2d Z_2 class. 
    """
    symm  = "Z_2"
    d     = 2
    alpha = 0
    beta  = 1/8
    gamma = 7/4
    delta = 15
    nu    = 1
    eta   = 1/4
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name


# 2d Potts q=3/ hard hexagon model 
class S3_2d(UniversalityClass):
    """ 
    Exact solution for 2d S_3 class from Baxter "Exactly Solved Models in Statistical Mechanics"
    """
    symm  = "S_3"
    d     = 2
    alpha = 1/3 
    beta  = 1/9
    gamma = 13/9 
    delta = 14
    nu    = 5/6
    eta   = 4/15
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name


# 2d Potts q=4/ Ashkin-Teller model
class S4_2d(UniversalityClass):
    """ 
    Exact solution for 2d S_4 class from Baxter, "Exactly Solved Models in Statistical Mechanics"
    """
    symm  = "S_4"
    d     = 2
    alpha = 2/3 
    beta  = 1/12
    gamma = 7/6
    delta = 15
    nu    = 2/3
    eta   = 1/4
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name


def reweight(X, pRW, p0, S):
    """ 
    Reweight an observable X computed at a simulation point p0 to a nearby
    simulation point pRW. We assume the action depends linearly on the simulation
    parameter, i.e. S' ~ p S

    By the way, if you are going to pair this with a jackknife to estimate e.g.
    where a response function is maximized, make sure you pick enough reweighting
    points to properly resolve where the maximum is.

    Args:
        X (np.array): Measurements to reweight. 
        pRW (float): Reweight to this target. 
        p0 (float): Simulation point.
        S (np.array): Measurements of the action (extensive) divided by parameter p. 
    """
    checkType(X,'array')
    checkType(S,'array')
    Z_i = np.exp( (pRW-p0)*(S-std_mean(S)) )
    Z   = np.sum(Z_i)
    return np.sum(X*Z_i/Z)

