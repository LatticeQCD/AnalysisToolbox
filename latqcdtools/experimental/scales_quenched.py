#
# scales_quenched.py
#
# L. Mazur, D. Clarke
#
# A collection of scales and related functions for pure SU(3) configurations.
#
import numpy as np
from latqcdtools.experimental.tools import fm_to_GeVinv, GeVinv_to_fm
from latqcdtools.experimental.scales_hisq import r1_MILC_2010
import latqcdtools.base.logger as logger



''' r0 scales '''
r0_phys_GeV= 0.469/0.1973269718  # TODO: add source


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
    r0_fm = 1.5092*r1_MILC_2010("fm") # about 0.469
    if units=="fm":
      return r0_fm
    elif units=="GeVinv":
      return fm_to_GeVinv(r0_fm)
    else:
      logger.TBError("Invalid unit specification for r0.")


def a_r0_invGeV(beta):
    return r0_hQCD_2014("GeVinv")/r0_div_a(beta)


def a_r0_fm(beta):
    return GeVinv_to_fm(a_r0_invGeV(beta))



''' t0 scales '''



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
