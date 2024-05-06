# 
# runningCoupling.py                                                               
# 
# D. Clarke
# 
# Python implementation of functions related to running coupling. Aims toward QCD, but implemented for
# SU(Nc) with arbitrary Nc and arbitrary Nf. The coefficients are given in terms of g(mu),
#
#     beta = -g^2 ( b0 g + b1 g^3 + ... ) 
#


import numpy as np
from scipy.special import zeta
from latqcdtools.base.check import checkType


ZETA_3 = zeta(3)


def CF(Nc) -> float:
    """ Casimir operator of defining representation of SU(Nc).

    Args:
        Nc (int): Number of colors.

    Returns:
        float : CF 
    """
    checkType(Nc,int)
    return (Nc**2-1)/(2*Nc)


def CA(Nc) -> float:
    """ Casimir operator of adjoint representation of SU(Nc).

    Args:
        Nc (int): Number of colors.

    Returns:
        float : CA 
    """
    checkType(Nc,int)
    return Nc 


# 10.1103/PhysRevLett.30.1346, 10.1103/PhysRevLett.30.1343
def b0(Nf,Nc=3) -> float:
    """ Universal one-loop beta coefficient for SU(Nc), expansion in g.  

    Args:
        Nf (int): Number of active quark flavors. 
        Nc (int): Number of colors.

    Returns:
        float : b0 
    """
    checkType(Nf,int)
    return ( 11*CA(Nc)/3 - 2*Nf/3 )


def b1(Nf,Nc=3) -> float:
    """ Universal two-loop beta coefficient for SU(Nc), expansion in g. 

    Args:
        Nf (int): Number of active quark flavors.
        Nc (int): Number of colors.

    Returns:
        float : b1 
    """
    checkType(Nf,int)
    return ( 34*CA(Nc)**2/3 - 2*CF(Nc)*Nf - 10*CA(Nc)*Nf/3 )


# Larin and Vermaseren, Phys Lett B 303 (1993) 334-336
def b2_dimreg_MSbar(Nf,Nc=3) -> float:
    """ Three-loop beta coefficient for SU(Nc) using dimensional regularization in MS-bar scheme. 

    Args:
        Nf (int): Number of active quark flavors.
        Nc (int): Number of colors.

    Returns:
        float : b2_MS-bar (dim reg) 
    """
    checkType(Nf,int)
    return ( 2857*CA(Nc)**3/54 + Nf   *(CF(Nc)**2-205*CF(Nc)*CA(Nc)/18-1415*CA(Nc)**2/54)
                               + Nf**2*(11*CF(Nc)/9+79*CA(Nc)/54)
            )


# van Ritbergen, Larin, and Vermaseren, Phys Lett B 400 (1997) 379-384
def b3_dimreg_MSbar(Nf,Nc=3) -> float:
    """ Four-loop beta coefficient for SU(Nc) using dimensional regularization in MS-bar scheme. 

    Args:
        Nf (int): Number of active quark flavors.
        Nc (int): Number of colors.

    Returns:
        float : b3_MS-bar (dim reg) 
    """
    checkType(Nf,int)
    return ( CA(Nc)**4*(150653/486-44*ZETA_3/9) + Nc**2*(Nc**2+36)/24*(704*ZETA_3/3-80/9)
                + Nf   *(CA(Nc)**3/2*(136*ZETA_3/3-39143/81) 
                          + CA(Nc)**2*CF(Nc)/2*(7073/243-656*ZETA_3/9) 
                          + CA(Nc)*CF(Nc)**2/2*(352*ZETA_3/9-4204/27) 
                          + 46*CF(Nc)**3/2 
                          + Nc*(Nc**2+6)/48*(512/9-1664*ZETA_3/3)
                         )
                + Nf**2*(CA(Nc)**2/4*(224*ZETA_3/9+7930/81)
                          + CF(Nc)**2/4*(1352/27-704*ZETA_3/9)
                          + CA(Nc)*CF(Nc)/4*(17152/243+448*ZETA_3/9)
                          + (Nc**4-6*Nc**2+18)/(96*Nc**2)*(512*ZETA_3/3-704/9)
                         )
                + Nf**3*(424*CA(Nc)/(8*243)+1232*CF(Nc)/(8*243)
                         )
            )


def beta_func(beta,Nf=3) -> float:
    """ QCD asymptotic scaling relation to two loops.

    Args:
        beta (float-like)
        Nf (int, optional): Number of active quark flavors. Defaults to 3.

    Returns:
        float: f_as(beta) 
    """
    checkType(Nf,int)
    B0 = b0(Nf)/(4*np.pi)**2
    B1 = b1(Nf)/(4*np.pi)**4
    return (B0*10/beta)**(-B1/(2*B0**2)) * np.exp(-beta/(20*B0))
