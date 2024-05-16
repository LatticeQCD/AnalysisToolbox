#
# referenceScales.py
#
# A collection of scales and related functions for pure SU(3) configurations.
#

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.math.polynomials import Polynomial, Rational
from latqcdtools.physics.runningCoupling import beta_func, b0, b1


CHECKBETARANGE = True


# Most recent update of scales in physical units
CY_phys  = { 
            'fk' : 2019,
            'r1' : 2010,
            'r0' : 2014,
            't0' : 2015 
            }


# Most recent update of scale parameterizations
CY_param = { 
            'fk' : 2021,
            'r1' : 2021,
            'r0' : 2017,
            't0' : 2015 
            }


def ignoreBetaRange():
    """ Turn off the beta range warnings. """
    global CHECKBETARANGE
    CHECKBETARANGE = False
    logger.warn('Squelching beta range warnings.')


def _betaRangeWarn(beta, beta_range):
    """ Many of these ansätze a(beta) have coefficients that were determined by performing a fit within a certain
    beta range. This warning flashes whenever you are using information outside of that range, where the ansatz
    is less likely be to reliable.

    Args:
        beta (float) or numpy array
        beta_range (array-like): min and max beta of range, in that order 
    """
    global CHECKBETARANGE
    print(beta)
    if  isinstance(beta , (float,np.floating)) :
        beta    = np.array([ beta , beta ]) 
    if CHECKBETARANGE:
        if np.sort(beta)[0] < beta_range[0] or np.sort(beta)[-1] > beta_range[1]:
            logger.warn("beta out of fit range [" + str(beta_range[0]) + "," + str(beta_range[1]) + "]",frame=3)


# -------------------------------------------------------------------------------------------------------- FITTING FORMS


def fit_2014Eos_eqB2(beta,c0,c2,d2):
    return Rational([0,c0,0,c2*10/beta],[1,0,d2*10/beta])(beta_func(beta))


def fit_tayloraLambda(beta,a,b,c):
    return Polynomial([0,a,0,b,0,c])(beta_func(beta))


# ----------------------------------------------------------------------------------------------- 2+1 FLAVOR HISQ SCALES


def allton_type_ansatz(beta, c0, c2, d2):
    return (c0 * beta_func(beta) + c2 * (10 / beta) * beta_func(beta) ** 3) / \
           (1 + d2 * (10 / beta) * beta_func(beta) ** 2)


# ===================================================== f_K scales


def a_times_fk(beta: float, year):
    """ Get a*f_k(beta).

    Args:
        beta (float)
        year (int): year that parameterization was determined 

    Returns:
        float: a*f_k 
    """
    # https://arxiv.org/pdf/2107.10011.pdf, 10.1103/PhysRevD.104.074512
    if str(year) == "2021":
        beta_range = [6.175, 7.220]
        _betaRangeWarn(beta, beta_range)
        c0fk = 7.486
        c2fk = 41935.0
        d2fk = 3273.0

    # 10.1103/PhysRevD.100.094510
    elif str(year) == "2014":
        beta_range = [6., 7.373]
        _betaRangeWarn(beta, beta_range)
        c0fk = 7.49415
        c2fk = 46049.0
        d2fk = 3671.0

    # https://arxiv.org/pdf/1111.1710 , 10.1103/PhysRevD.85.054503 
    elif str(year) == "2012":
        beta_range = [6., 6.8]
        _betaRangeWarn(beta, beta_range)
        c0fk = 7.65667
        c2fk = 32911.0
        d2fk = 2388.0

    else:
        logger.TBError("No fit parameters for ", str(year))

    return allton_type_ansatz(beta, c0fk, c2fk, d2fk)


# ====================================================== r1 scales


def a_div_r1(beta, year):
    """ Get a/r_1(beta).

    Args:
        beta (float)
        year (int): year that parameterization was determined 

    Returns:
        float: a/r_1
    """
    # https://arxiv.org/pdf/2107.10011.pdf, 10.1103/PhysRevD.104.074512
    if str(year) == "2021":
        c0 = 43.16
        c2 = 339472
        d2 = 5452.0
    # https://arxiv.org/pdf/1710.05024.pdf
    elif str(year) == "2018":
        beta_range = [7.030, 8.4]
        _betaRangeWarn(beta, beta_range)
        c0 = 43.12
        c2 = 347008
        d2 = 5584
    # https://arxiv.org/pdf/1407.6387.pdf
    elif str(year) == "2014":
        c0 = 43.1
        c2 = 343236.0
        d2 = 5514.0
    # https://arxiv.org/pdf/1111.1710.pdf
    elif str(year) == "2012":
        c0 = 44.06
        c2 = 272102.0
        d2 = 4281.0
    else:
        logger.TBError("No fit parameters for year", str(year))
    return allton_type_ansatz(beta, c0, c2, d2)


# ================================================= strange quark mass: line of constant physics


# fit taken from 1407.6387v2
def r1_times_ms_2014(beta):
    nf   = 3
    mRGI = 0.2609
    m1   = 35600
    m2   = -21760
    m3   = 2.67*10**7
    dm1  = 2420
    num  = 1 +  m1*(10.0/beta)*beta_func(beta)**2+m2*(10.0/beta)**2*beta_func(beta)**2 + m3*(10.0/beta)*beta_func(beta)**4
    den  = 1 + dm1*(10.0/beta)*beta_func(beta)**2
    return (20*b0(nf)/beta)**(4/9)*mRGI*num/den


def a_times_ms_2014(beta):
    return r1_times_ms_2014(beta) * a_div_r1(beta, "2014")


# ------------------------------------------------------------------------------------------------------ QUENCHED SCALES


def wuppertal_type_ansatz(beta,c1,c2,c3,c4):
    Nf = 0
    B0 = b0(Nf)/(4*np.pi)**2
    B1 = b1(Nf)/(4*np.pi)**4
    return np.exp( (beta/(12*B0)+B1/(2.*B0**2)*np.log(6*B0/beta))
                   * (1+c1/beta+c2/beta**2)/(1+c3/beta+c4/beta**2) )


# ===================================================== r0 scales


def r0_div_a(beta,year):
    """ Get r0/a(beta).

    Args:
        beta (float)
        year (int): year that parameterization was determined 

    Returns:
        float: r0/a 
    """
    # https://arxiv.org/abs/1709.07612
    if str(year) == "2017":
        c1 = -8.9664
        c2 = 19.21
        c3 = -5.25217
        c4 = 0.606828
    # https://arxiv.org/abs/1503.05652
    elif str(year) == "2015":
        beta_range = [5.7, 6.92]
        _betaRangeWarn(beta, beta_range)
        c1 = -8.17273
        c2 = 14.9600
        c3 = -3.95983
        c4 = -5.30334
    else:
        logger.TBError("No fit parameters for year", str(year))
    return wuppertal_type_ansatz(beta,c1,c2,c3,c4)


# ===================================================== t0 scales


# Based on https://arxiv.org/pdf/1503.05652.pdf
# Latest update at 2017/01/11 by Lukas Mazur
def sqrtt0_div_a(beta):
    """ Get sqrt(t0/a)(beta)

    Args:
        beta (float)

    Returns:
        float: sqrt(t0/a) 
    """
    c1 = -9.945
    c2 = 24.191
    c3 = -5.334
    c4 = 1.452
    return wuppertal_type_ansatz(beta,c1,c2,c3,c4)
