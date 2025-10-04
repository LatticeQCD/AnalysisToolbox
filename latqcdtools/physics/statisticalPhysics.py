# 
# statisticalPhysics.py                                                               
# 
# D. Clarke
# 
# A collection of methods relevant for statistical physics calculations. This includes
# tables for critical exponents and temperatures from the literature. We assume that
# estimates for multiple critical parameters coming from the same paper have correctly
# accounted for correlations; i.e. we do not here account for the fact that e.g. nu and
# eta coming from the same Markov chain are correlated.
#


import numpy as np
import gvar as gv
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.statistics.statistics import gaudif 


def _printExponent(prefix, exponent):
    if exponent is not None:
        logger.info(prefix, exponent)


def _getParameter(parameterInfo):
    """ 
    Give the weighted average of literature values of some parameter. These are 
    in units with k_B=J=1 and no magnetic field. If there is only one known value, 
    give that back. 
    """
    if len(parameterInfo)==1:
        for key in parameterInfo:
            return parameterInfo[key]
    else:
        _params = []
        for key in parameterInfo:
            _params.append(parameterInfo[key])
        _params = np.array(_params)
        return np.mean(_params) 


def _statisticalConsistencyCheck(parameterInfo):
    for key1 in parameterInfo:
        for key2 in parameterInfo:
            m1 = gv.mean(parameterInfo[key1])
            e1 = gv.sdev(parameterInfo[key1])
            m2 = gv.mean(parameterInfo[key2])
            e2 = gv.sdev(parameterInfo[key2])
            if gaudif(m1,e1,m2,e2)<0.05:
                logger.warn("Statistical tension between",key1,key2)


def _compareWithZero(err,tol) -> bool:
    return ( np.abs(gv.mean(err))<=gv.sdev(err) ) or (err<tol)



class UniversalityClass:

    """ 
    Skeleton universality class from which all others inherit.
    """

    symm  = None
    d     = None

    alphas = {} 
    betas  = {} 
    gammas = {} 
    deltas = {} 
    nus    = {} 
    etas   = {} 
    omegas = {} 
    Tcs    = {}

    alpha = None
    beta  = None
    gamma = None
    delta = None
    nu    = None
    eta   = None
    omega = None
    Tc    = None

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

    def hyperscalingCheck(self, tol=1e-12) -> bool:
        lpass = True
        err1 = 2*self.beta+self.gamma-2+self.alpha
        err2 = 2*self.beta*self.delta-self.gamma-2+self.alpha
        err3 = self.nu*self.d-2+self.alpha
        lpass *= _compareWithZero(err1,tol) 
        lpass *= _compareWithZero(err2,tol) 
        lpass *= _compareWithZero(err3,tol) 
        if not lpass:
            logger.TBFail(f"{self.name()} fails at least one hypersclaing relation")
            logger.TBFail(f"err1 = {err1}, err2 = {err2}, err3 = {err3}")
        return lpass 


# 3d XY model
class O2_3d(UniversalityClass):
    symm = "O(2)"
    d    = 3
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    etas = {
        '10.1007/JHEP08(2016)036' : gv.gvar('0.03852(64)'),
    }
    nus  = {
        '10.1007/JHEP08(2016)036' : gv.gvar('0.6719(11)'),
    }

    eta   = _getParameter(etas)
    nu    = _getParameter(nus)
    gamma = nu*(2-eta)
    beta  = 0.5*(nu*d-gamma)
    alpha = 2 - nu*d
    delta = (gamma+2*alpha)/(2*beta)


# 3d Heisenberg model
class O3_3d(UniversalityClass):
    symm = "O(3)"
    d    = 3
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    etas = {
        '10.1007/JHEP08(2016)036'  : gv.gvar('0.0386(12)'),
    }
    nus  = {
        '10.1007/JHEP08(2016)036'  : gv.gvar('0.7121(28)'),
    }
    Tcs  = {
        '10.1103/PhysRevB.43.6087' : gv.gvar('1.44321(21)'),
        '10.1103/PhysRevB.48.936'  : gv.gvar('1.44300(21)')
    }
    _statisticalConsistencyCheck(Tcs)

    Tc    = _getParameter(Tcs)
    eta   = _getParameter(etas)
    nu    = _getParameter(nus)
    gamma = nu*(2-eta)
    beta  = 0.5*(nu*d-gamma)
    alpha = 2 - nu*d
    delta = (gamma+2*alpha)/(2*beta)


class O4_3d(UniversalityClass):
    symm = "O(4)"
    d    = 3
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    nus    = {
        '10.1103/PhysRevD.51.2404'         : gv.gvar('0.7479(90)')
    }
    betas  = {
        '10.1103/PhysRevD.51.2404'         : gv.gvar('0.3836(46)')
    }
    gammas = {
        '10.1103/PhysRevD.51.2404'         : gv.gvar('1.477(18)')
    }
    deltas = {
        '10.1016/j.nuclphysb.2003.09.060'  : gv.gvar('4.824(9)'),
    }

    beta  = _getParameter(betas) 
    delta = _getParameter(deltas) 
    gamma = _getParameter(gammas) 
    nu    = _getParameter(nus) 
    alpha = 2.-beta*(1.+delta)


# 3d Ising model
class Z2_3d(UniversalityClass):
    symm = "Z_2"
    d    = 3
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    etas = {
        '10.1007/s10955-014-1042-7'    : gv.gvar('0.03631(3)'),
        '10.1007/JHEP08(2016)036'      : gv.gvar('0.0362978(20)'),
    }
    nus  = {
        '10.1007/s10955-014-1042-7'    : gv.gvar('0.62999(5)'),
        '10.1007/JHEP08(2016)036'      : gv.gvar('0.629971(4)'),
    }
    Tcs  = {
        '10.1103/PhysRevB.29.4030'     : gv.gvar('4.51154(13)'),
        '10.1103/PhysRevB.32.1720'     : gv.gvar('4.51162(11)'),
        '10.1016/0378-4371(94)90490-1' : gv.gvar('4.511463(45)'),
        '10.1103/PhysRevB.82.174433'   : gv.gvar('4.5115232(17)')
    }
    _statisticalConsistencyCheck(etas)
    _statisticalConsistencyCheck(nus)
    _statisticalConsistencyCheck(Tcs)

    eta   = _getParameter(etas)
    nu    = _getParameter(nus)
    alpha = 2 - nu*d
    gamma = nu*(2 - eta)
    beta  = (nu*d -gamma)/2
    delta = nu*d/beta-1.


# 2d Ising model
class Z2_2d(UniversalityClass):
    """ 
    Onsager solution for 2d Z_2 class. 
    """
    symm = "Z_2"
    d    = 2
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    alpha = 0
    beta  = 1/8
    gamma = 7/4
    delta = 15
    nu    = 1
    eta   = 1/4
    Tc    = 2/np.log(1+np.sqrt(2))


# 2d Potts q=3/ hard hexagon model 
class S3_2d(UniversalityClass):
    """ 
    Exact solution for 2d S_3 class from Baxter "Exactly Solved Models in Statistical Mechanics"
    """
    symm = "S_3"
    d    = 2
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    alpha = 1/3 
    beta  = 1/9
    gamma = 13/9 
    delta = 14
    nu    = 5/6
    eta   = 4/15


# 2d Potts q=4/ Ashkin-Teller model
class S4_2d(UniversalityClass):
    """ 
    Exact solution for 2d S_4 class from Baxter, "Exactly Solved Models in Statistical Mechanics"
    """
    symm = "S_4"
    d    = 2
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    alpha = 2/3 
    beta  = 1/12
    gamma = 7/6
    delta = 15
    nu    = 2/3
    eta   = 1/4


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

