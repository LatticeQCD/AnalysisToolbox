# 
# constants.py
# 
# D. Clarke 
# 
# Many constants and unit conversions that are relevant for physics.
#

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType


# Base constants for unit conversions
hceVm           = 197.3269788e-9   # hbar c in [eV m]. PDG 2018. DOI: 10.1103/PhysRevD.98.030001.
cms             = 299792458        # c in [m/s]. NIST 2018.
days_per_year   = 365.2422         # NIST 2018.
days_per_month  = days_per_year/12
meters_per_mile = 1609.344         # NIST 2023 based on international foot (not survey foot).
BTU_per_Wh      = 3.412            # Many possible definitions, don't take too seriously. A stupid unit, indeed.
kBJdivK         = 1.380649e-23     # kB in [J/K]. NIST 2018.
eC              = 1.602176634e-19  # e in [C]. NIST 2018.


# List of scientific prefixes
_prefix = { "Q"  : 1e30,
            "R"  : 1e27,
            "Y"  : 1e24,  
            "Z"  : 1e21, 
            "E"  : 1e18, 
            "P"  : 1e15, 
            "T"  : 1e12, 
            "G"  : 1e9, 
            "M"  : 1e6, 
            "k"  : 1e3, 
            "h"  : 1e2, 
            "da" : 1e1, 
            "1"  : 1,
            1    : 1,
            "d"  : 1e-1,
            "c"  : 1e-2, 
            "m"  : 1e-3, 
            "u"  : 1e-6, 
            "n"  : 1e-9, 
            "p"  : 1e-12, 
            "f"  : 1e-15, 
            "a"  : 1e-18, 
            "z"  : 1e-21, 
            "y"  : 1e-24, 
            "r"  : 1e-27, 
            "q"  : 1e-30  }


# Alphabetical order for easier finding
_baseUnits = [
             "BTU",    # British thermal unit
             "degC",   # Celsius
             "degF",   # Fahrenheit
             "eV",     # electron-volt
             "ft",     # feet
             "h",      # hour
             "J",      # Joule
             "K",      # Kelvin
             "m",      # meter
             "m/s",    # meter/second
             "mi",     # mile
             "mi/h",   # mile/hour
             "min",    # minute
             "s",      # second
             "d",      # day
             "W",      # Watt
             "Wh",     # Watt-hour
             "Wh/y",   # Watt-hour/year
             "y",      # year
             ]


def _separatePrefix(units):
    checkType(units,str)
    if units in _baseUnits: 
        prefix=1
        baseUnit=units
    else:
        prefix=units[0]
        baseUnit=units[1:]
    return prefix, baseUnit


def convert(x,unit1,unit2) -> float:
    """ General method for doing unit conversions.

    Args:
        x (float): measurement in [unit1]. 
        unit1 (str): Original units.
        unit2 (str): Target units.

    Returns:
        float: measurement in [unit2]. 
    """
    p1, u1 = _separatePrefix(unit1)
    p2, u2 = _separatePrefix(unit2)
    if not p1 in _prefix:
        logger.TBError('Unknown prefix',p1)
    if not p2 in _prefix:
        logger.TBError('Unknown prefix',p2)
    u1u2 = (u1,u2)

    if u1.endswith('inv'):
        num=1/_prefix[p1]
    else:
        num=_prefix[p1]
    if u2.endswith('inv'):
        den=1/_prefix[p2]
    else:
        den=_prefix[p2]
    fac = num/den

    if u1==u2:
        result = x

    # time 
    elif u1u2==('s','min'):
        result = x/60
    elif u1u2==('min','s'):
        result = x*60
    elif u1u2==('s','h'):
        result = x/60**2
    elif u1u2==('h','s'):
        result = x*60**2
    elif u1u2==('s','y'):
        result = x/(days_per_year*24*60**2)
    elif u1u2==('y','s'):
        result = x*(days_per_year*24*60**2)
    elif u1u2==('h','y'):
        result = x/(days_per_year*24)
    elif u1u2==('y','h'):
        result = x*(days_per_year*24)
    elif u1u2==('d','y'):
        result = x/days_per_year
    elif u1u2==('y','d'):
        result = x*days_per_year

    # distance
    elif u1u2==('mi','m'):
        result = x*meters_per_mile
    elif u1u2==('m','mi'):
        result = x/meters_per_mile
    elif u1u2==('mi','ft'):
        result = x*5280
    elif u1u2==('ft','mi'):
        result = x/5280
    elif u1u2==('ft','m'):
        result = convert( convert(x,'ft','mi'), 'mi','m')
    elif u1u2==('m','ft'):
        result = convert( convert(x,'m','mi'), 'mi','ft')
        
    # speed
    elif u1u2==('mi/h','m/s'):
        result = convert( convert(x,'mi','m'), 's','h')
    elif u1u2==('m/s','mi/h'):
        result = convert( convert(x,'m','mi'), 'h','s')

    # power 
    elif u1u2==('W','Wh/y'):
        result = convert(x,'y','h')
    elif u1u2==('Wh/y','W'):
        result = convert(x,'h','y')

    # temperature
    elif u1u2==('K','degC'):
        result = x - 273.15/fac 
    elif u1u2==('degC','K'):
        result = x + 273.15/fac 
    elif u1u2==('degC','degF'):
        result = x*9/5 + 32/fac
    elif u1u2==('degF','degC'):
        result = (5/9)*(x-32/fac) 
    elif u1u2==('degF','K'):
        result = convert( convert(x,'degF','degC'), 'degC','K' )
    elif u1u2==('K','degF'):
        result = convert( convert(x,'K','degC'), 'degC','degF' )

    # energy
    elif u1u2==('Wh','BTU'):
        result = x*BTU_per_Wh
    elif u1u2==('BTU','Wh'):
        result = x/BTU_per_Wh
    elif u1u2==('Wh','J'):
        result = convert(x,'h','s')
    elif u1u2==('J','Wh'):
        result = convert(x,'s','h')
    elif u1u2==('eV','J'): # [J] = [V C]
        result = x*eC 
    elif u1u2==('J','eV'):
        result = x/eC 

    # natural units
    elif u1u2==('m','eVinv'):
        result = x/hceVm
    elif u1u2==('eVinv','m'):
        result = x*hceVm
    elif u1u2==('eV','minv'):
        result = x/hceVm
    elif u1u2==('minv','eV'):
        result = x*hceVm
    elif u1u2==('eV','K'):
        result = convert(x,'eV','J')/kBJdivK
    elif u1u2==('K','eV'):
        result = convert(x*kBJdivK,'J','eV')

    else:
        logger.TBError('No rule for conversion of ['+u1+'] to ['+u2+']') 

    return fac*result


hcMeVfm = convert( convert(hceVm,"eV","MeV"), "m","fm" )
hcGeVfm = convert(hcMeVfm,"MeV","GeV")


def fm_to_MeVinv(x) -> float:
    return x/hcMeVfm
def fm_to_GeVinv(x) -> float:
    return x/hcGeVfm

def MeV_to_fminv(x) -> float:
    return x/hcMeVfm
def GeV_to_fminv(x) -> float:
    return x/hcGeVfm

def MeVinv_to_fm(x) -> float:
    return hcMeVfm*x
def GeVinv_to_fm(x) -> float:
    return hcGeVfm*x


# ------------------------------------------------------------------------------------------------------ PARTICLE MASSES 


def M_mu_phys(year=2020,units="MeV",returnErr=False):
    """ Physical value of the muon mass. """
    if year==2020:
        # PDG 2020. DOI: https://doi.org/10.1093/ptep/ptaa104.
        m_mu_MeV, m_mu_MeV_err = 105.6583745, 0.0000024
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return convert(m_mu_MeV,"MeV",units), convert(m_mu_MeV_err,"MeV",units)
    else:
        return convert(m_mu_MeV,"MeV",units)


def M_pi0_phys(year=2022,units="MeV",returnErr=False):
    """ Physical value of the pi0 mass. """
    if year==2022:
        # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097
        m_pi0_MeV, m_pi0_MeV_err = 134.9768, 0.0005
    elif year==2014:
        # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
        m_pi0_MeV, m_pi0_MeV_err = 134.9766, 0.0006
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return convert(m_pi0_MeV,"MeV",units), convert(m_pi0_MeV_err,"MeV",units)
    else:
        return convert(m_pi0_MeV,"MeV",units)


def M_pipm_phys(year=2022,units="MeV",returnErr=False):
    """ Physical value of the pi+/- mass. """
    if year==2022:
        # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097
        m_pipm_MeV, m_pipm_MeV_err = 139.57039, 0.00018
    elif year==2014:
        # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
        m_pipm_MeV, m_pipm_MeV_err = 139.57018, 0.00035
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return convert(m_pipm_MeV,"MeV",units), convert(m_pipm_MeV_err,"MeV",units)
    else:
        return convert(m_pipm_MeV,"MeV",units)


def M_rho_phys(year=2022,units="MeV",returnErr=False):
    """ Physical value of the rho mass. """
    if year==2022:
        # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097
        m_rho_MeV, m_rho_MeV_err = 775.26, 0.23
    elif year==2014:
        # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
        m_rho_MeV, m_rho_MeV_err = 775.26, 0.25
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return convert(m_rho_MeV,"MeV",units), convert(m_rho_MeV_err,"MeV",units)
    else:
        return convert(m_rho_MeV,"MeV",units)


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
    return convert( fkMeV/np.sqrt(2.), "MeV", units )


def frho_phys(year=2017,units="GeV",returnErr=False):
    """ Physical value of the rho decay constant. """
    if year==2017:
        # HPQCD 2017. DOI: https://doi.org/10.1103/PhysRevD.93.014503. Figure 6.
        frhoGeV, frhoGeV_err = 0.21, 0.01
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return convert( frhoGeV, "GeV", units ), convert( frhoGeV_err, "GeV", units )
    else:
        return convert( frhoGeV, "GeV", units )
    

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
        return convert(LMS,"MeV",units), convert(LMSerr,"MeV",units)
    else:
        return convert(LMS,"MeV",units)


def Rproton_phys(year=2018,units="fm",returnErr=False):
    """ Physical value of proton charge radius. """
    if year==2018:
        # From NIST 2018. 
        R, Rerr =  0.8414, 0.0019
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return convert(R,"fm",units), convert(Rerr,"fm",units)
    else:
        return convert(R,"fm",units)


def r1_phys(year=2010,units="fm",returnErr=False):
    """ Physical value of Sommer scale r1. """    
    if year==2010:
        # r1 taken from MILC 2010. arXiv:1012.0868. 
        r1fm, r1fm_err = 0.3106, np.sqrt( 0.0008**2 + 0.0014**2 + 0.0004**2 )
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return convert(r1fm,"fm",units), convert(r1fm_err,"fm",units)
    else:
        return convert(r1fm,"fm",units)


def r0_phys(year=2014,units="fm",returnErr=False):
    """ Physical value of Sommer scale r0. """    
    if year==2014:
        # Use r0/r1 from hotQCD 2014. DOI: https://doi.org/10.1103/PhysRevD.90.094503.
        r1, r1_err = r1_phys(year=2010,units=units,returnErr=True)
        r0, r0_err = 1.5092*r1, 1.5092*r1_err
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return r0, r0_err 
    else:
        return r0 


def sqrt_t0_phys(year=2015,units="fm",returnErr=False):
    """ Gradient flow scale t0^(1/2). """
    if year==2015:
        # TODO: Where did this come from? Is the year correct? 
        sqrtt0_div_r0 = 0.334
        sqrtt0        = sqrtt0_div_r0 * r0_phys(year=2014,units=units,returnErr=False) 
    else:
        logger.TBError("Invalid year specification.")
    if returnErr:
        return sqrtt0, None 
    else:
        return sqrtt0

