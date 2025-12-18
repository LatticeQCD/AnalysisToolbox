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


CATCHTENSION = True


def ignoreTension():
    """ 
    Turn off statistical tension warnings. 
    """
    global CATCHTENSION
    CATCHTENSION = False
    logger.warn("Ignoring tests for statistical tension.")


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


def _statisticalConsistencyCheck(parameterInfo,obs):
    global CATCHTENSION
    if CATCHTENSION:
        for key1 in parameterInfo:
            for key2 in parameterInfo:
                m1 = gv.mean(parameterInfo[key1])
                e1 = gv.sdev(parameterInfo[key1])
                m2 = gv.mean(parameterInfo[key2])
                e2 = gv.sdev(parameterInfo[key2])
                if gaudif(m1,e1,m2,e2)<0.05:
                    logger.warn(f"Statistical tension for {obs} between",key1,key2)


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
        _printExponent("   eta =",self.eta)
        logger.info()

    def hyperscalingCheck(self, tol=1e-12) -> bool:
        lpass = True
        err1 = 2*self.beta+self.gamma-2+self.alpha
        err2 = 2*self.beta*self.delta-self.gamma-2+self.alpha
        err3 = self.nu*self.d-2+self.alpha
        err4 = self.gamma-self.nu*(2-self.eta)
        lpass *= _compareWithZero(err1,tol) 
        lpass *= _compareWithZero(err2,tol) 
        lpass *= _compareWithZero(err3,tol) 
        lpass *= _compareWithZero(err4,tol) 
        if not lpass:
            logger.TBFail(f"{self.name()} fails at least one hyperscaling relation")
            logger.TBFail(f"err1 = {err1}, err2 = {err2}, err3 = {err3}, err4 = {err4}")
        return lpass 


# 3d XY model
class O2_3d(UniversalityClass):
    symm = "O(2)"
    d    = 3
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    etas = {
        '10.1007/JHEP08(2016)036'     : gv.gvar('0.03852(64)'),
        '10.1103/PhysRevB.100.224517' : gv.gvar('0.03810(8)'),
    }
    nus  = {
        '10.1007/JHEP08(2016)036'     : gv.gvar('0.6719(11)'),
        '10.1103/PhysRevB.100.224517' : gv.gvar('0.67169(7)'),
    }

    eta   = _getParameter(etas)
    nu    = _getParameter(nus)
    gamma = nu*(2-eta)
    alpha = 2 - nu*d
    beta  = 0.5*(2-alpha-gamma)
    delta = (2-alpha+gamma)/(2*beta) 


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
    _statisticalConsistencyCheck(Tcs,'Tc')

    Tc    = _getParameter(Tcs)
    eta   = _getParameter(etas)
    nu    = _getParameter(nus)
    gamma = nu*(2-eta)
    alpha = 2 - nu*d
    beta  = 0.5*(2-alpha-gamma)
    delta = (2-alpha+gamma)/(2*beta) 


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
    eta   = 2.-gamma/nu


# O(N), N=infinity limit
class Oinf_3d(UniversalityClass):
    symm = "O(inf)"
    d    = 3
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name
    beta  = 1/2
    delta = 5
    alpha = 2.-beta*(1.+delta)
    nu    = (2.-alpha)/d
    gamma = 2. - alpha - 2*beta 
    eta   = 2. - gamma/nu


# 3d Ising model. There is a nice review article here: 10.1016/S0370-1573(00)00127-7.
# Some of the results were taken from the MCMC/FSS section of tables 2 and 3 of that
# reference. There they list Kc, yt, and yh. I used those values to compute Tc, nu,
# and beta, so there may be minor rounding issues. beta was extracted from yh taking
# correlations with nu into account.
#
# The results from 10.1103/PhysRevB.44.5081 are in statistical tension with many of
# the other results, and they cause a failure of the hyperscaling check, so they have 
# been excluded from averages.
class Z2_3d(UniversalityClass):
    symm = "Z_2"
    d    = 3
    def __repr__(self) -> str:
        return super().__repr__()+':'+self.name

    etas = {
        '10.1007/s10955-014-1042-7'    : gv.gvar('0.03631(3)'),      # 2014
        '10.1007/JHEP08(2016)036'      : gv.gvar('0.0362978(20)'),   # 2016
    }
    nus  = {
        '10.1016/0378-4371(94)90490-1' : gv.gvar('0.62893(79)'),     # 1994
        '10.1088/0305-4470/28/22/007'  : gv.gvar('0.63012(79)'),     # 1995
        '10.1142/S0129183196000247'    : gv.gvar('0.6309(28)'),      # 1996
        '10.1088/0305-4470/30/1/006'   : gv.gvar('0.6309(12)'),      # 1997
        '10.1103/PhysRevB.59.11471'    : gv.gvar('0.62980(52)'),     # 1999
        '10.1088/0305-4470/32/26/304'  : gv.gvar('0.62960(20)'),     # 1999
        '10.1142/S0129183199000929'    : gv.gvar('0.63032(56)'),     # 1999
        '10.1088/0305-4470/32/1/004'   : gv.gvar('0.62941(48)'),     # 1999
        '10.1007/s10955-014-1042-7'    : gv.gvar('0.62999(5)'),      # 2014
        '10.1007/JHEP08(2016)036'      : gv.gvar('0.629971(4)'),     # 2016
    }
    betas = {
        '10.1016/0378-4371(94)90490-1' : gv.gvar('0.3258(44)'),      # 1994
        '10.1088/0305-4470/28/22/007'  : gv.gvar('0.3267(10)'),      # 1995
        '10.1142/S0129183196000247'    : gv.gvar('0.3237(24)'),      # 1996
        '10.1103/PhysRevB.59.11471'    : gv.gvar('0.32643(37)'),     # 1999
        '10.1088/0305-4470/32/26/304'  : gv.gvar('0.32607(16)'),     # 1999
        '10.1142/S0129183199000929'    : gv.gvar('0.32682(38)'),     # 1999
        '10.1088/0305-4470/32/1/004'   : gv.gvar('0.32597(40)'),     # 1999
    }
    Tcs  = {
        '10.1103/PhysRevB.29.4030'     : gv.gvar('4.51154(13)'),     # 1984
        '10.1103/PhysRevB.32.1720'     : gv.gvar('4.51162(11)'),     # 1985
        '10.1209/0295-5075/16/2/003'   : gv.gvar('4.511528(20)'),    # 1991
        '10.1143/JPSJ.60.1978'         : gv.gvar('4.511475(61)'),    # 1991
        '10.1016/0378-4371(93)90208-L' : gv.gvar('4.511658(81)'),    # 1993
        '10.1016/0378-4371(94)90490-1' : gv.gvar('4.511463(45)'),    # 1994
        '10.1088/0305-4470/28/22/007'  : gv.gvar('4.511524(20)'),    # 1995
        '10.1088/0305-4470/29/17/042'  : gv.gvar('4.511528(12)'),    # 1996
        '10.1142/S0129183196000247'    : gv.gvar('4.511516(20)'),    # 1996
        '10.1103/PhysRevE.54.2291'     : gv.gvar('4.51152(31)'),     # 1996
        '10.1142/S0129183199000929'    : gv.gvar('4.511524(20)'),    # 1999
        '10.1103/PhysRevB.82.174433'   : gv.gvar('4.5115232(17)'),   # 2010  
    }
    _statisticalConsistencyCheck(etas,'eta')
    _statisticalConsistencyCheck(nus,'nu')
    _statisticalConsistencyCheck(betas,'beta')
    _statisticalConsistencyCheck(Tcs,'Tc')

    Tc    = _getParameter(Tcs)
    eta   = _getParameter(etas)
    nu    = _getParameter(nus)
    beta  = _getParameter(betas)
    alpha = 2 - nu*d
    gamma = nu*(2 - eta)
    delta = nu*d/beta-1.
    eta   = 2.-gamma/nu


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

