# 
# constants.py
# 
# D. Clarke 
# 
# Many constants that are relevant for physics
#

import numpy as np
import latqcdtools.base.logger as logger


# Taken from PDG 2018. DOI: 10.1103/PhysRevD.98.030001.
hcMeVfm=197.3269788
hcGeVfm=hcMeVfm/1000


def fm_to_MeVinv(x):
    return x/hcMeVfm
def fm_to_GeVinv(x):
    return x/hcGeVfm

def MeV_to_fminv(x):
    return x/hcMeVfm
def GeV_to_fminv(x):
    return x/hcGeVfm

def MeVinv_to_fm(x):
    return hcMeVfm*x
def GeVinv_to_fm(x):
    return hcGeVfm*x


def gethc(units="MeVfm"):
    if units=="MeVfm":
        return hcMeVfm
    elif units=="GeVfm":
        return hcGeVfm
    else:
        logger.TBError("Invalid unit specification.")


def MeVtoUnits(value,name,units):
    if units == "MeV":
        return value
    elif units == "GeV":
        return value/1000.
    elif units == "fminv":
        return MeV_to_fminv(value)
    else:
        logger.TBError("Invalid unit specification for " + name + ".")


def GeVtoUnits(value,name,units):
    if units == "GeV":
        return value
    elif units == "MeV":
        return value*1000.
    elif units == "fminv":
        return GeV_to_fminv(value)
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


# ------------------------------------------------------------------------------------------------------ PARTICLE MASSES 


def M_mu_phys(year=2020,units="MeV",returnErr=False):
    """ Physical value of the muon mass. """
    name = "M_mu"
    if year==2020:
        # PDG 2020. DOI: https://doi.org/10.1093/ptep/ptaa104.
        m_mu_MeV, m_mu_MeV_err = 105.6583745, 0.0000024
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return MeVtoUnits(m_mu_MeV,name,units), MeVtoUnits(m_mu_MeV_err,name,units)
    else:
        return MeVtoUnits(m_mu_MeV,name,units)


def M_pi0_phys(year=2022,units="MeV",returnErr=False):
    """ Physical value of the pi0 mass. """
    name = "M_pi^0"
    if year==2022:
        # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097
        m_pi0_MeV, m_pi0_MeV_err = 134.9768, 0.0005
    elif year==2014:
        # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
        m_pi0_MeV, m_pi0_MeV_err = 134.9766, 0.0006
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return MeVtoUnits(m_pi0_MeV,name,units), MeVtoUnits(m_pi0_MeV_err,name,units)
    else:
        return MeVtoUnits(m_pi0_MeV,name,units)


def M_pipm_phys(year=2022,units="MeV",returnErr=False):
    """ Physical value of the pi+/- mass. """
    name = "M_pi^+-"
    if year==2022:
        # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097
        m_pipm_MeV, m_pipm_MeV_err = 139.57039, 0.00018
    elif year==2014:
        # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
        m_pipm_MeV, m_pipm_MeV_err = 139.57018, 0.00035
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return MeVtoUnits(m_pipm_MeV,name,units), MeVtoUnits(m_pipm_MeV_err,name,units)
    else:
        return MeVtoUnits(m_pipm_MeV,name,units)


def M_rho_phys(year=2022,units="MeV",returnErr=False):
    """ Physical value of the rho mass. """
    name = "M_rho"
    if year==2022:
        # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097
        m_rho_MeV, m_rho_MeV_err = 775.26, 0.23
    elif year==2014:
        # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
        m_rho_MeV, m_rho_MeV_err = 775.26, 0.25
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return MeVtoUnits(m_rho_MeV,name,units), MeVtoUnits(m_rho_MeV_err,name,units)
    else:
        return MeVtoUnits(m_rho_MeV,name,units)


# ------------------------------------------------------------------------------------------------------ DECAY CONSTANTS


def fk_phys(year=2019,units="MeV"):
    """ Physical value of Kaon decay constant. """
    if year==2019:
        # FLAG 2019. DOI: 10.1140/epjc/s10052-019-7354-7. Section 4.6.
        fkMeV = 155.7
    elif year==2018:
        # PDG 2018. DOI: 10.1103/PhysRevD.98.030001. Section 84.5.1.
        fkMeV = 155.72
    elif year==2012:
        # PDG 2012. DOI: 10.1103/PhysRevD.86.010001. Page 949 under meson listings.
        fkMeV = 156.1
    else:
        logger.TBError("Invalid year specification.")
    return MeVtoUnits( fkMeV/np.sqrt(2.), "fK", units )


def frho_phys(year=2017,units="GeV",returnErr=False):
    """ Physical value of the rho decay constant. """
    if year==2017:
        # HPQCD 2017. DOI: https://doi.org/10.1103/PhysRevD.93.014503. Figure 6.
        frhoGeV, frhoGeV_err = 0.21, 0.01
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return GeVtoUnits( frhoGeV, "frho", units ), GeVtoUnits( frhoGeV_err, "frhoerr", units )
    else:
        return GeVtoUnits( frhoGeV, "frho", units )
    

# ------------------------------------------------------------------------------------------------------ OTHER CONSTANTS 


def alpha_e(year=2018,returnErr=False):
    """ Fine structure constant. """
    # NIST 2018 CODATA recommended value.
    if year==2018:
        alpha, alpha_err  = 7.2973525693e-3, 0.0000000011e-3
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return alpha, alpha_err 
    else:
        return alpha 


def lambda_MSbar_phys(year=2021,units="MeV",returnErr=False):
    """ Physical value of MS-bar lambda parameter. """
    if year==2021:
        # Kaon decay constant taken from FLAG 2021. arXiv: 2111.09849
        LMS, LMSerr = 339, 12
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return MeVtoUnits(LMS,"lambda-MSbar",units), MeVtoUnits(LMSerr,"lambda-MSbarerr",units)
    else:
        return MeVtoUnits(LMS,"lambda-MSbar",units)