#
# referenceScales.py
#
# A collection of scales and related functions for pure SU(3) configurations.
#

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.physics.unitConversions import fm_to_GeVinv, GeVinv_to_fm, fm_to_MeVinv, MeV_to_fminv
from latqcdtools.math.polynomials import Polynomial,Rational


def beta_func(beta):
    """ Asymptotic scaling relation for Nf=3 up to two loops. """
    nf = 3
    b0 = (11 - 2 * nf / 3) / (4 * np.pi) ** 2
    b1 = (102 - 38 * nf / 3) / (4 * np.pi) ** 4
    return (b0 * 10 / beta) ** (-b1 / (2 * b0 ** 2)) * np.exp(-beta / (20 * b0))


def MeVtoUnits(value,name,units):
    if units == "MeV":
        return value
    elif units == "fminv":
        return MeV_to_fminv(value)
    else:
        logger.TBError("Invalid unit specification for " + name + ".")


def fmtoUnits(value,name,units):
    if units == "fm":
        return value
    elif units == "MeVinv":
        return fm_to_MeVinv(value)
    elif units == "GeVinv":
        return fm_to_GeVinv(value)
    else:
        logger.TBError("Invalid unit specification for " + name + ".")


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


def print_out_of_beta_range_warning(beta, beta_range):
    if beta < beta_range[0] or beta > beta_range[1]:
        logger.warn("beta out of fit range [" + str(beta_range[0]) + "," + str(beta_range[1]) + "]")


def a_times_fk(beta: float, year, suppress_warnings=False):

    # https://arxiv.org/pdf/2107.10011.pdf, 10.1103/PhysRevD.104.074512
    if str(year) == "2021":
        beta_range = [6.175, 7.220]
        if not suppress_warnings:
            print_out_of_beta_range_warning(beta, beta_range)
        c0fk = 7.486
        c2fk = 41935.0
        d2fk = 3273.0

    # 10.1103/PhysRevD.100.094510
    elif str(year) == "2014":
        # TODO add beta range
        c0fk = 7.49415
        c2fk = 46049.0
        d2fk = 3671.0

    # TODO add source
    elif str(year) == "2012":
        # TODO add beta range
        c0fk = 7.65667
        c2fk = 32911.0
        d2fk = 2388.0

    else:
        logger.TBError("No fit parameters for ", str(year))

    return allton_type_ansatz(beta, c0fk, c2fk, d2fk)


def a_fk_invGeV(beta: float, year, suppress_warnings=False):
    if str(year) == "2021":
        fKexp = fk_phys(2019,"MeV")
    elif str(year) == "2014" or str(year) == "2012":
        fKexp = fk_phys(2012,"MeV")
    else:
        logger.TBError("No fit parameters for ", str(year))
    return (a_times_fk(beta, year, suppress_warnings) * 1000) / fKexp


def a_fk_fm(beta, year):
    return GeVinv_to_fm(a_fk_invGeV(beta, year))


def fk_phys(year=2019,units="MeV"):
    """ Physical value of Kaon decay constant. """
    if year==2019:
        # Kaon decay constant taken from FLAG 2019. DOI: 10.1140/epjc/s10052-019-7354-7. Section 4.6.
        fkMeV = 155.7 / np.sqrt(2.)
    elif year==2018:
        # Kaon decay constant taken from PDG 2018. DOI: 10.1103/PhysRevD.98.030001. Section 84.5.1.
        fkMeV = 155.72 / np.sqrt(2.)
    elif year==2012:
        # Kaon decay constant taken from PDG 2012. DOI: 10.1103/PhysRevD.86.010001. Page 949 under meson listings.
        fkMeV = 156.1 / np.sqrt(2.)
    else:
        logger.TBError("Invalid year specification for physical value of fK.")
    return MeVtoUnits(fkMeV,"fK",units)



# ====================================================== r1 scales


def a_div_r1(beta, year, suppress_warnings=False):

    # https://arxiv.org/pdf/2107.10011.pdf, 10.1103/PhysRevD.104.074512
    if str(year) == "2021":
        # TODO add beta range
        c0 = 43.16
        c2 = 339472
        d2 = 5452.0

    # https://arxiv.org/pdf/1710.05024.pdf
    elif str(year) == "2018":
        beta_range = [7.030, 8.4]
        if not suppress_warnings:
            print_out_of_beta_range_warning(beta, beta_range)
        c0 = 43.1
        c2 = 343236.0
        d2 = 5514.0

    # https://arxiv.org/pdf/1407.6387.pdf
    elif str(year) == "2014":
        # TODO add beta range
        c0 = 43.1
        c2 = 343236.0
        d2 = 5514.0

    # https://arxiv.org/pdf/1111.1710.pdf
    elif str(year) == "2012":
        # TODO add beta range
        c0 = 44.06
        c2 = 272102.0
        d2 = 4281.0

    else:
        logger.TBError("No fit parameters for year", str(year))
    return allton_type_ansatz(beta, c0, c2, d2)


def a_r1_invGeV(beta, year, suppress_warnings=False):
    return fm_to_GeVinv(r1_MILC_2010("fm") * a_div_r1(beta, year, suppress_warnings))


def a_r1_fm(beta, year, suppress_warnings=False):
    return r1_MILC_2010("fm") * a_div_r1(beta, year, suppress_warnings)


def r1_MILC_2010(units):
    """ r1 taken from MILC 2010. arXiv:1012.0868. """
    r1fm = 0.3106
    return fmtoUnits(r1fm,"r1",units)


def r1err_MILC_2010(units):
    """ Error bar from MILC 2010. arXiv:101.0868. """
    r1errfm = np.sqrt( 0.0008**2 + 0.0014**2 + 0.0004**2 )
    return fmtoUnits(r1errfm,"r1_err",units)


# ================================================= strange quark mass: line of constant physics


# fit taken from 1407.6387v2
def r1_times_ms_2014(beta):
    nf = 3
    b0 = (11 - 2 * nf / 3) / (4 * np.pi) ** 2
    mRGI = 0.2609
    m1 = 35600
    m2 = -21760
    m3 = 2.67 * 10 ** 7
    dm1 = 2420
    num = 1 + m1 * (10.0 / beta) * beta_func(beta) ** 2 + m2 * (10.0 / beta) ** 2 * beta_func(beta) ** 2 + m3 * (10.0 / beta) * beta_func(beta) ** 4
    den = 1 + dm1 * (10.0 / beta) * beta_func(beta) ** 2
    return (20 * b0 / beta) ** (4.0 / 9.0) * mRGI * num / den


def a_times_ms_2014(beta):
    return r1_times_ms_2014(beta) * a_div_r1(beta, "2014")


# ------------------------------------------------------------------------------------------------------ QUENCHED SCALES


# ===================================================== r0 scales


# Fit ansatz from https://arxiv.org/pdf/1503.05652.pdf
# Coefficients from https://arxiv.org/abs/1709.07612
def r0_div_a(beta):
    b0=11./(4*np.pi)**2
    b1=102/(4*np.pi)**4
    c1=-8.9664
    c2=19.21
    c3=-5.25217
    c4=0.606828
    return np.exp( (beta/(12*b0)+b1/(2.*b0**2)*np.log(6*b0/beta))
                   * (1+c1/beta+c2/beta**2)/(1+c3/beta+c4/beta**2) )


# Use r0/r1 from hotQCD 2014. DOI: https://doi.org/10.1103/PhysRevD.90.094503.
# Use r1=0.3106 from MILC. arXiv:1012.0868.
def r0_hQCD_2014(units):
    r0fm = 1.5092*r1_MILC_2010("fm") # about 0.469
    return fmtoUnits(r0fm,"r0",units)


def r0err_hQCD_2014(units):
    r0errfm = 1.5092*np.sqrt( 0.0008**2 + 0.0014**2 + 0.0004**2 )
    return fmtoUnits(r0errfm,"r0err",units)


def a_r0_invGeV(beta):
    return r0_hQCD_2014("GeVinv")/r0_div_a(beta)


def a_r0_fm(beta):
    return GeVinv_to_fm(a_r0_invGeV(beta))


# ===================================================== t0 scales


# Based on https://arxiv.org/pdf/1503.05652.pdf
# Latest update at 2017/01/11 by Lukas Mazur
def sqrtt0(beta):
    b0=11./(4*np.pi)**2
    b1=102/(4*np.pi)**4
    c1=-9.945
    c2=24.191
    c3=-5.334
    c4=1.452
    return np.exp( (beta/(12*b0)+b1/(2.*b0**2)*np.log(6*b0/beta))
                   * (1+c1/beta+c2/beta**2)/(1+c3/beta+c4/beta**2) )


sqrtt0r0_cont = 0.334
sqrtt0_phys   = sqrtt0r0_cont*r0_hQCD_2014("GeVinv")


def a_t0_invGeV(beta):
    return sqrtt0_phys/sqrtt0(beta)


def a_t0_fm(beta):
    return GeVinv_to_fm(a_t0_invGeV(beta))


# --------------------------------------------------------------------------------------------------------- OTHER SCALES


def lambda_MSbar_phys(year=2021,units="MeV",returnErr=False):
    """ Physical value of MS-bar lambda parameter. """
    if year==2021:
        # Kaon decay constant taken from FLAG 2021. arXiv: 2111.09849
        LMS, LMSerr = 339, 12
    else:
        logger.TBError("Invalid year specification for physical value of lambda-MSbar.")
    if returnErr:
        return MeVtoUnits(LMS,"lambda-MSbar",units), MeVtoUnits(LMSerr,"lambda-MSbarerr",units)
    else:
        return MeVtoUnits(LMS,"lambda-MSbar",units)


# ------------------------------------------------------------------------------------------------------ LEGACY WRAPPERS


def fk_FLAG_2019(units):
    return fk_phys(2019,units)
def fk_PDG_2018(units):
    return fk_phys(2018,units)
def fk_PDG_2012(units):
    return fk_phys(2012,units)
