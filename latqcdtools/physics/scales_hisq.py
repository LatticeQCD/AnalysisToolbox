#
# scales_hisq.py
#
# L. Mazur, D. Clarke
#
# A collection of scales and related functions for Nf=2+1 HISQ configurations.
#
import numpy as np
from latqcdtools.physics.unitConversions import GeVinv_to_fm, MeV_to_fminv, fm_to_GeVinv, fm_to_MeVinv
import latqcdtools.base.logger as logger


def beta_func(beta):
    nf = 3
    b0 = (11 - 2 * nf / 3) / (4 * np.pi) ** 2
    b1 = (102 - 38 * nf / 3) / (4 * np.pi) ** 4
    return (b0 * 10 / beta) ** (-b1 / (2 * b0 ** 2)) * np.exp(-beta / (20 * b0))


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
    # TODO add source
    if str(year) == "2021":
        fKexpnew = 155.7
    # TODO add source
    elif str(year) == "2014" or str(year) == "2012":
        fKexpnew = 156.1
    else:
        logger.TBError("No fit parameters for ", str(year))
    return (a_times_fk(beta, year, suppress_warnings) * np.sqrt(2.) * 1000) / fKexpnew


def a_fk_fm(beta, year):
    return GeVinv_to_fm(a_fk_invGeV(beta, year))


# Experimental Kaon decay constant taken from PDG 2018. DOI: 10.1103/PhysRevD.98.030001
# It's in section 84.5.1.
def fk_PDG_2018(units):
    fkMeV = 155.72 / np.sqrt(2.)
    if units == "MeV":
        return fkMeV
    elif units == "fminv":
        return MeV_to_fminv(fkMeV)
    else:
        logger.TBError("Invalid unit specification for fk.")


# Experimental Kaon decay constant taken from PDG 2012. DOI: 10.1103/PhysRevD.86.010001
# It's on page 949 under meson particle listings.
def fk_PDG_2012(units):
    fkMeV = 156.1 / np.sqrt(2.)
    if units == "MeV":
        return fkMeV
    elif units == "fminv":
        return MeV_to_fminv(fkMeV)
    else:
        logger.TBError("Invalid unit specification for fk.")


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
        logger.TBError("No fit parameters for ", str(year))
    return allton_type_ansatz(beta, c0, c2, d2)


def a_r1_invGeV(beta, year, suppress_warnings=False):
    return fm_to_GeVinv(r1_MILC_2010("fm") * a_div_r1(beta, year, suppress_warnings))


def a_r1_fm(beta, year, suppress_warnings=False):
    return r1_MILC_2010("fm") * a_div_r1(beta, year, suppress_warnings)


# r1 taken from MILC 2010. arXiv:1012.0868.
def r1_MILC_2010(units):
    r1fm = 0.3106
    if units == "fm":
        return r1fm
    elif units == "MeVinv":
        return fm_to_MeVinv(r1fm)
    elif units == "GeVinv":
        return fm_to_GeVinv(r1fm)
    else:
        logger.TBError("Invalid unit specification for r1.")


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
