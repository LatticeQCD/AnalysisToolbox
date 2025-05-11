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
from latqcdtools.statistics.statistics import std_mean, weighted_mean, weighted_variance
from latqcdtools.base.printErrorBars import getValuesFromErrStr
from latqcdtools.base.utilities import toNumpy


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
    Tcs   = {} # J=1, no magnetic field

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

    def Tc(self,paper=None):
        """ Give the weighted average of literature values of Tc. These are in units with
        k_B=J=1 and no magnetic field. If there is only one known value, give that back.
        The 2d, Z(2) univerality class just gives back the Onsager solution. If you want
        to access a particular result, you can pass the DOI of the paper. """ 
        if len(self.Tcs)==0:
            logger.TBRaise('No Tc data for this universality class.')
        elif paper is not None:
            return getValuesFromErrStr(self.Tcs[paper])
        elif (self.symm=="Z_2") and (self.d==2):
            return self.Tcs['Onsager']
        elif len(self.Tcs)==1:
            for key in self.Tcs:
                return getValuesFromErrStr(self.Tcs[key])
        else:
            _Tcs, _Tces = [],[]
            for key in self.Tcs:
                m, e = getValuesFromErrStr(self.Tcs[key])
                _Tcs.append(m)
                _Tces.append(e)
            _Tcs, _Tces = toNumpy(_Tcs,_Tces)
            return weighted_mean(_Tcs, _Tces), np.sqrt(weighted_variance(_Tces))


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
    Tcs = {
        '10.1103/PhysRevB.43.6087' : '1.44321(21)',
        '10.1103/PhysRevB.48.936'  : '1.44300(21)'
    }
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
    Tcs = {
        '10.1103/PhysRevB.29.4030'     : '4.51154(13)',
        '10.1103/PhysRevB.32.1720'     : '4.51162(11)',
        '10.1016/0378-4371(94)90490-1' : '4.511463(45)',
        '10.1103/PhysRevB.82.174433'   : '4.5115232(17)'
    }
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
    Tcs = {
        'Onsager' : 2/np.log(1+np.sqrt(2))
    }
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
    checkType(np.ndarray,X=X)
    checkType(np.ndarray,S=S)
    arg = (pRW-p0)*S
    # Subtracting by np.max(arg) is a trick to avoid overflow and improve
    # numerical stability.
    Z_i = np.exp( arg-np.max(arg) ) 
    Z   = np.sum(Z_i)
    return np.sum(X*Z_i/Z)

