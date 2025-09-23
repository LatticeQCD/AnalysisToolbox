# 
# constants.py
# 
# D. Clarke 
# 
# Many constants and unit conversions that are relevant for physics. These values do not necessarily
# represent the most up-to-date or appropriate choices; rather they are a collection of constants I
# have needed for work or teaching classes. This is especially true of constants like particle masses,
# where choosing one value over another may be subtle. 
#

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.math.math import quadrature


# Base constants for unit conversions
hceVm           = 197.3269788e-9   # hbar c in [eV m]. PDG 2018. DOI: 10.1103/PhysRevD.98.030001.
cms             = 299792458        # c in [m/s]. NIST 2018.
days_per_year   = 365.2422         # NIST 2018.
days_per_month  = days_per_year/12
meters_per_mile = 1609.344         # NIST 2023 based on international foot (not survey foot).
BTU_per_Wh      = 3.412            # Many possible definitions--don't take too seriously. A stupid unit, indeed.
kBJdivK         = 1.380649e-23     # kB in [J/K]. NIST 2018.
eC              = 1.602176634e-19  # e in [C]. NIST 2018.


# See corresponding _phys() functions for the references
_fkerrs2012   = np.array([0.2   ,0.8   ,0.2   ])
_fkerrs2018   = np.array([0.17  ,0.45  ,0.16  ])
_fpierrs2018  = np.array([0.01  ,0.03  ,0.13  ])
_r1fmerrs2010 = np.array([0.0008,0.0014,0.0004])


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
    checkType(str,units=units)
    if units in _baseUnits: 
        prefix=1
        baseUnit=units
    else:
        prefix=units[0]
        baseUnit=units[1:]
    return prefix, baseUnit


def convert(x,unit1,unit2) -> float:
    """ 
    General method for doing unit conversions. He knows about scientific prefixes like G, M, and so on.
    If the unit ends in 'inv', it is interpreted as 1/unit.

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
        logger.TBRaise('Unknown prefix',p1)
    if not p2 in _prefix:
        logger.TBRaise('Unknown prefix',p2)
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
        logger.TBRaise('No rule for conversion of ['+u1+'] to ['+u2+']') 

    return fac*result


# This block is to support some legacy code.
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
def fminv_to_MeV(x) -> float:
    return x*hcMeVfm


class physicalConstant():

    name = None

    def __init__(self,name,scale,units):
        """
        Wrap a dictionary of scale values in physical units.

        Args:
            name (str)
            scale (dict): A dictionary of lists indexed by year. 
            units (str): Base units. See convert method for what units are allowed. 
        """
        checkType(str,name=name)
        checkType(dict,scale=scale)
        if units is not None:
            checkType(str,units=units)
        for key in scale:
            checkType(tuple,scale=scale[key])
        self.name=name
        self.scale=scale
        self.scaleUnits=units
    
    def __repr__(self) -> str:
        if self.name is None:
            return 'physicalConstant'
        return self.name 

    def getValue(self,year,units,returnErr,normalize=1.):
        """
        Retrieve the value of this scale.

        Args:
            year (int): Year this was measured. 
            units (str): Get it in these units. Use None for unitless quantities. 
            returnErr (bool): Return both mean and uncertainty as tuple.
            normalize (real): Return scale/normalize. Defaults to 1.
        """
        checkType('int',year=year)
        if units is not None:
            checkType(str,units=units)
        checkType(bool,returnErr=returnErr)
        checkType("scalar",normalize=normalize)
        if not year in self.scale:
            logger.TBRaise(f"Invalid year specification {year}. Allowed years:",list(self.scale.keys()))
        if self.scaleUnits is None:
            if units is not None:
                logger.TBRaise("Tried to convert a unitless quantity to units",units)
            val = self.scale[year][0]
            err = self.scale[year][1]
        else:
            val = convert(self.scale[year][0],self.scaleUnits,units)/normalize
            err = convert(self.scale[year][1],self.scaleUnits,units)/normalize
        if returnErr:
            return val, err 
        else:
            return val 


# ------------------------------------------------------------------------------------------------------ PARTICLE MASSES 


def M_u_phys(year=2024,units="MeV",returnErr=False):
    scale = {
        2024: (2.16, 0.07), # PDG 2024. DOI: 10.1103/PhysRevD.110.030001
    }
    return physicalConstant("M_u",scale,"MeV").getValue(year,units,returnErr) 


def M_d_phys(year=2024,units="MeV",returnErr=False):
    scale = {
        2024: (4.7, 0.07), # PDG 2024. DOI: 10.1103/PhysRevD.110.030001
    }
    return physicalConstant("M_d",scale,"MeV").getValue(year,units,returnErr) 


def M_s_phys(year=2024,units="MeV",returnErr=False):
    scale = {
        2024: (93.5, 0.8), # PDG 2024. DOI: 10.1103/PhysRevD.110.030001
    }
    return physicalConstant("M_s",scale,"MeV").getValue(year,units,returnErr) 


def M_c_phys(year=2024,units="GeV",returnErr=False):
    scale = {
        2024: (1.273,0.0046), # PDG 2024. DOI: 10.1103/PhysRevD.110.030001
    }
    return physicalConstant("M_c",scale,"MeV").getValue(year,units,returnErr) 


def M_b_phys(year=2024,units="GeV",returnErr=False):
    scale = {
        2024: (4.183,0.007), # PDG 2024. DOI: 10.1103/PhysRevD.110.030001
    }
    return physicalConstant("M_b",scale,"MeV").getValue(year,units,returnErr) 


def M_t_phys(year=2024,units="GeV",returnErr=False):
    scale = {
        2024: (172.57,0.29), # PDG 2024. DOI: 10.1103/PhysRevD.110.030001. direct measurement
    }
    return physicalConstant("M_t",scale,"MeV").getValue(year,units,returnErr) 


def M_mu_phys(year=2022,units="MeV",returnErr=False):
    scale = {
        2022: (105.6583755, 0.0000023), # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097.
        2020: (105.6583745, 0.0000024)  # PDG 2020. DOI: https://doi.org/10.1093/ptep/ptaa104.
    }
    return physicalConstant("M_mu",scale,"MeV").getValue(year,units,returnErr) 


def M_pi0_phys(year=2022,units="MeV",returnErr=False):
    scale = {
        2022: (134.9768, 0.0005), # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097.
        2014: (134.9766, 0.0006)  # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
    }
    return physicalConstant("M_pi0",scale,"MeV").getValue(year,units,returnErr) 


def M_pipm_phys(year=2022,units="MeV",returnErr=False):
    scale = {
        2022: (139.57039, 0.00018), # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097.
        2014: (139.57018, 0.00035)  # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
    }
    return physicalConstant("M_pi+/-",scale,"MeV").getValue(year,units,returnErr) 


def M_K0_phys(year=2022,units="MeV",returnErr=False):
    scale = {
        2022: (497.611, 0.013), # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097.
    }
    return physicalConstant("M_K0",scale,"MeV").getValue(year,units,returnErr) 


def M_Kpm_phys(year=2022,units="MeV",returnErr=False):
    scale = {
        2022: (493.677, 0.013), # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097.
    }
    return physicalConstant("M_K+/-",scale,"MeV").getValue(year,units,returnErr) 


def M_rho_phys(year=2022,units="MeV",returnErr=False):
    scale = {
        2022: (775.26, 0.23), # PDG 2022. DOI: https://doi.org/10.1093/ptep/ptac097.
        2014: (775.26, 0.25)  # PDG 2014. DOI: 10.1088/1674-1137/38/9/090001
    }
    return physicalConstant("M_rho",scale,"MeV").getValue(year,units,returnErr) 


def M_proton_phys(year=2024,units="MeV",returnErr=False):
    scale = {
        2024: (938.27208816, 0.00000029), # PDG 2024. DOI: 10.1103/PhysRevD.110.030001
    }
    return physicalConstant("M_proton",scale,"MeV").getValue(year,units,returnErr) 


def M_neutron_phys(year=2024,units="MeV",returnErr=False):
    scale = {
        2024: (939.56542052, 0.00000054), # PDG 2024. DOI: 10.1103/PhysRevD.110.030001
    }
    return physicalConstant("M_neutron",scale,"MeV").getValue(year,units,returnErr) 



# ------------------------------------------------------------------------------------------------------ DECAY CONSTANTS



def fk_phys(year=2019,units="MeV",returnErr=False):
    """ 
    Physical value of Kaon decay constant, f_K+/-. Scaled by sqrt(2.). 
    """
    scale = {
        2019: (155.7 , 0.7),                     # FLAG 2019. DOI: 10.1140/epjc/s10052-019-7354-7. Section 4.6, Nf=2+1.
        2018: (155.72, quadrature(_fkerrs2018)), # PDG 2018. DOI: 10.1103/PhysRevD.98.030001. Section 84.5.1.
        2012: (156.1 , quadrature(_fkerrs2012)), # PDG 2012. DOI: 10.1103/PhysRevD.86.010001. Page 949 under meson listings.
    }
    return physicalConstant("f_K+/-",scale,"MeV").getValue(year,units,returnErr,normalize=np.sqrt(2.)) 


def fpi_phys(year=2018,units="MeV",returnErr=False):
    """
    Physical value of the pion decay constant, f_pi+/-. 
    """
    scale = {
        2018: (130.50, quadrature(_fpierrs2018)), # PDG 2018. DOI: 10.1103/PhysRevD.98.030001. Section 84.5.1.
    }
    return physicalConstant("f_pi+/-",scale,"MeV").getValue(year,units,returnErr) 


def frho_phys(year=2017,units="GeV",returnErr=False):
    """ 
    Physical value of the rho decay constant. 
    """
    scale = {
        2017: (0.21, 0.01), # HPQCD 2017. DOI: https://doi.org/10.1103/PhysRevD.93.014503. Figure 6.
    }
    return physicalConstant("f_rho",scale,"GeV").getValue(year,units,returnErr) 


# ---------------------------------------------------------------------------------------------------- LATTICE CONSTANTS 


def w0_phys(year=2013,units="fm",returnErr=False):
    """ 
    Gradient flow scale w0.
    """    
    scale = {
        2013: (0.1715, 0.0009), # w0 taken from HPQCD. DOI: 10.1103/PhysRevD.88.074504. Eq (18).
    }
    return physicalConstant("w0",scale,"fm").getValue(year,units,returnErr) 


def r1_phys(year=2010,units="fm",returnErr=False):
    """ 
    Physical value of Sommer scale r1. 
    """    
    scale = {
        2010: (0.3106, quadrature(_r1fmerrs2010)), # r1 taken from MILC 2010. arXiv:1012.0868. 
    }
    return physicalConstant("r1",scale,"fm").getValue(year,units,returnErr) 


_r1phys2010, _r1phys2010e = r1_phys(year=2010,units="fm",returnErr=True)
def r0_phys(year=2014,units="fm",returnErr=False):
    """ 
    Physical value of Sommer scale r0. 
    """    
    scale = {
        2014: (_r1phys2010*1.5092, _r1phys2010e*1.5092), # Use r0/r1 from hotQCD 2014. DOI: https://doi.org/10.1103/PhysRevD.90.094503.
    }
    return physicalConstant("r0",scale,"fm").getValue(year,units,returnErr) 


# ------------------------------------------------------------------------------------------------------ OTHER CONSTANTS 


def alpha_e(year=2018,returnErr=False):
    """ 
    Fine structure constant. 
    """
    scale = {
        2018: (7.2973525693e-3, 0.0000000011e-3), # NIST 2018 CODATA recommended value.
    }
    return physicalConstant("alpha_e",scale,None).getValue(year,None,returnErr) 


def lambda_MSbar_phys(year=2021,units="MeV",returnErr=False):
    """ 
    Physical value of MS-bar lambda parameter. 
    """
    scale = {
        2021: ( 339, 12 ), # Taken from FLAG 2021. arXiv: 2111.09849
    }
    return physicalConstant("lambda_MSbar",scale,"MeV").getValue(year,units,returnErr) 


def Rproton_phys(year=2018,units="fm",returnErr=False):
    """ 
    Physical value of proton charge radius. 
    """
    scale = {
        2018: ( 0.8414, 0.0019 ), # Taken from FLAG 2021. arXiv: 2111.09849
    }
    return physicalConstant("Rproton",scale,"fm").getValue(year,units,returnErr) 
