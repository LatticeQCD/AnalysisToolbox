#
# referenceScales.py
#
# A collection of scales and related functions for pure SU(3) configurations.
#

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.physics.constants import convert, fm_to_GeVinv, GeVinv_to_fm, fk_phys
from latqcdtools.math.polynomials import Polynomial, Rational
from latqcdtools.physics.betaFunction import beta_func, b0, b1


CHECKBETARANGE = True


# Current Year for various scales and parameterizations
CY_A_TIMES_FK = 2021
CY_FK_PHYS    = 2019
CY_A_DIV_R1   = 2021
CY_A_DIV_R0   = 2017


def ignoreBetaRange():
    """ Turn off the beta range warnings. """
    global CHECKBETARANGE
    CHECKBETARANGE = False
    logger.warn('Squelching beta range warnings.')


def print_out_of_beta_range_warning(beta, beta_range):
    """ Many of these ansätze a(beta) have coefficients that were determined by performing a fit within a certain
    beta range. This warning flashes whenever you are using information outside of that range, where the ansatz
    is less likely be to reliable.

    Args:
        beta (float)
        beta_range (array-like): min and max beta of range, in that order 
    """
    global CHECKBETARANGE
    if CHECKBETARANGE:
        if beta < beta_range[0] or beta > beta_range[1]:
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
        print_out_of_beta_range_warning(beta, beta_range)
        c0fk = 7.486
        c2fk = 41935.0
        d2fk = 3273.0

    # 10.1103/PhysRevD.100.094510
    elif str(year) == "2014":
        c0fk = 7.49415
        c2fk = 46049.0
        d2fk = 3671.0

    # TODO add source
    elif str(year) == "2012":
        c0fk = 7.65667
        c2fk = 32911.0
        d2fk = 2388.0

    else:
        logger.TBError("No fit parameters for ", str(year))

    return allton_type_ansatz(beta, c0fk, c2fk, d2fk)


def a_fk_invGeV(beta: float, year):
    if str(year) == "2021":
        fKexp = fk_phys(2019,"MeV")
    elif str(year) == "2014" or str(year) == "2012":
        fKexp = fk_phys(2012,"MeV")
    else:
        logger.TBError("No fit parameters for ", str(year))
    return (a_times_fk(beta, year) * 1000) / fKexp


def a_fk_fm(beta, year):
    return GeVinv_to_fm(a_fk_invGeV(beta, year))


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
        print_out_of_beta_range_warning(beta, beta_range)
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


def a_r1_invGeV(beta, year):
    return fm_to_GeVinv(r1_MILC_2010("fm") * a_div_r1(beta, year))


def a_r1_fm(beta, year):
    return r1_MILC_2010("fm") * a_div_r1(beta, year)


def r1_MILC_2010(units):
    """ r1 taken from MILC 2010. arXiv:1012.0868. """
    r1fm = 0.3106
    return convert(r1fm,"fm",units)


def r1err_MILC_2010(units):
    """ Error bar from MILC 2010. arXiv:101.0868. """
    r1errfm = np.sqrt( 0.0008**2 + 0.0014**2 + 0.0004**2 )
    return convert(r1errfm,"fm",units)


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
    nf = 0
    return np.exp( (beta/(12*b0(nf))+b1(nf)/(2.*b0(nf)**2)*np.log(6*b0(nf)/beta))
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
        c1=-8.9664
        c2=19.21
        c3=-5.25217
        c4=0.606828
    # https://arxiv.org/abs/1503.05652
    elif str(year) == "2015":
        beta_range = [5.7, 6.92]
        print_out_of_beta_range_warning(beta, beta_range)
        c1=-8.17273
        c2=14.9600
        c3=-3.95983
        c4=-5.30334
    else:
        logger.TBError("No fit parameters for year", str(year))
    return wuppertal_type_ansatz(beta,c1,c2,c3,c4)


# Use r0/r1 from hotQCD 2014. DOI: https://doi.org/10.1103/PhysRevD.90.094503.
# Use r1=0.3106 from MILC. arXiv:1012.0868.
def r0_hQCD_2014(units):
    r0fm = 1.5092*r1_MILC_2010("fm") # about 0.469
    return convert(r0fm,"fm",units)


def r0err_hQCD_2014(units):
    r0errfm = 1.5092*np.sqrt( 0.0008**2 + 0.0014**2 + 0.0004**2 )
    return convert(r0errfm,"fm",units)


def a_r0_invGeV(beta):
    return r0_hQCD_2014("GeVinv")/r0_div_a(beta)


def a_r0_fm(beta):
    return GeVinv_to_fm(a_r0_invGeV(beta))


# ===================================================== t0 scales


# Based on https://arxiv.org/pdf/1503.05652.pdf
# Latest update at 2017/01/11 by Lukas Mazur
def sqrtt0(beta):
    nf = 0
    c1=-9.945
    c2=24.191
    c3=-5.334
    c4=1.452
    return np.exp( (beta/(12*b0(nf))+b1(nf)/(2.*b0(nf)**2)*np.log(6*b0(nf)/beta))
                   * (1+c1/beta+c2/beta**2)/(1+c3/beta+c4/beta**2) )


sqrtt0r0_cont = 0.334
sqrtt0_phys   = sqrtt0r0_cont*r0_hQCD_2014("GeVinv")


def a_t0_invGeV(beta):
    return sqrtt0_phys/sqrtt0(beta)


def a_t0_fm(beta):
    return GeVinv_to_fm(a_t0_invGeV(beta))
