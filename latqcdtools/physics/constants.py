# 
# constants.py
# 
# D. Clarke 
# 
# Many constants that are relevant for physics
#


import latqcdtools.base.logger as logger


# Taken from PDG 2018. DOI: 10.1103/PhysRevD.98.030001.
hcMeVfm=197.3269788
hcGeVfm=hcMeVfm/1000


def fm_to_MeVinv(x):
    return x/hcMeVfm
def fm_to_GeVinv(x):
    return x/hcGeVfm
def MeVinv_to_fm(x):
    return hcMeVfm*x
def GeVinv_to_fm(x):
    return hcGeVfm*x
def MeV_to_fminv(x):
    return x/hcMeVfm


def gethc(units="MeVfm"):
    if units=="MeVfm":
        return hcMeVfm
    elif units=="GeVfm":
        return hcGeVfm
    else:
        logger.TBError("Invalid unit specification.")


# Fine structure constant. NIST 2018 CODATA recommended value.
alpha      = 7.2973525693e-3
alpha_err  = 0.0000000011e-3


# PDG 2020.
m_mu_MeV     = 105.6583745
m_mu_MeV_err =   0.0000024