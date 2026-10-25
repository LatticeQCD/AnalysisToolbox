# 
# correlators.py                                                               
# 
# D. Clarke
# 
# Some tools useful when analyzing correlators 
# 
import numpy as np
from latqcdtools.base.check import checkType, InvalidValueError
from latqcdtools.math.optimize import persistentSolve
from latqcdtools.base.utilities import unvector
import latqcdtools.base.logger as logger


def foldCorrelator(corr) -> np.ndarray:
    """
    On a periodic lattice, information for t>Nt/2 is redundant. This method
    combines corr below and above that threshold to improve statistics.
    
    Args:
        corr (np.ndarray)

    Returns:
        np.ndarray: folded correlator 
    """
    checkType(np.ndarray,corr=corr)
    if corr.ndim != 1:
        logger.TBRaise('corr should be 1-d')
    Nt          = len(corr)
    half_Nt     = Nt // 2
    first_half  = corr[1:half_Nt]
    second_half = corr[Nt-1:half_Nt:-1]
    averaged    = (first_half + second_half) / 2.0
    result      = np.concatenate(([corr[0]], averaged, [corr[half_Nt]]))
    return result


def _periodicLHS(corr,m,it,Nt):
    """
    Helper function for finding effective mass curve, baking in the fact
    that the lattice is periodic. See e.g. Gattringer and Lang eq. (6.57).
    Meant to be used with a solver.
    """
    rat1 = corr[it]/corr[it+1]
    rat2 = np.cosh(m*(it-Nt//2))/np.cosh(m*(it+1-Nt//2))
    return rat1-rat2


def effectiveMass(corr,algorithm='simple',guess=None,Nt=None) -> np.ndarray:
    """
    Get m_eff(t). Where this plateaus can be used as a first estimate for
    the ground state. See e.g. Gattringer and Lang eq. (6.56). The 'simple'
    algorithm ignores the backward-propagating part of the correlator, so it
    is biased low when exp(-m*(Nt-2t)) is not negligible; 'periodic' accounts for it.

    Args:
        corr (np.ndarray)
        algorithm (str): 'simple' or 'periodic'
        guess (float): initial guess for periodic mass meff[0]
        Nt (int, optional): temporal extent of the lattice. Needed with 'periodic' when corr
          is folded, i.e. has Nt//2+1 entries. Defaults to len(corr), i.e. an unfolded corr.

    Returns:
        np.ndarray: m_eff(t)
    """
    checkType(np.ndarray,corr=corr)
    checkType(str,algorithm=algorithm)
    if guess is not None:
        checkType('real',guess=guess)
        if algorithm=='simple':
            logger.TBRaise('simple algorithm requires no guess')
    else:
        if algorithm=='periodic':
            logger.TBRaise('periodic algorithm requires a starting guess for meff[0]')
    if Nt is None:
        Nt = len(corr)
    else:
        checkType('int',Nt=Nt)
        if algorithm=='simple':
            logger.TBRaise('simple algorithm does not use Nt')
        if len(corr) not in [Nt, Nt//2+1]:
            logger.TBRaise(f'len(corr)={len(corr)} should be Nt={Nt} or Nt//2+1 (folded)')
    meff = []
    if algorithm=='simple':
        try:
            for it in range(len(corr)-1):
                meff.append( np.log(corr[it]/corr[it+1]) )
        except InvalidValueError:
            logger.TBRaise('corr=',corr,exception=InvalidValueError)
    elif algorithm=='periodic':
        mit = guess
        for it in range(len(corr)-1):
            def LHS(m):
                return _periodicLHS(corr,m,it,Nt)
            # Depending on which solver converged, mit may be a scalar, 0-d array, or 1-element array
            mit = unvector(persistentSolve(LHS,guess=mit,maxiter=3000))
            meff.append( mit )
    else:
        logger.TBRaise('algorithm must be simple or periodic')
    return np.array(meff)
