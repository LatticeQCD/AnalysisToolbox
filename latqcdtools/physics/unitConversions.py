# 
# unitConversions.py                                                               
# 
# D. Clarke 
# 
# Common unit conversions in lattice QCD. 
#

import latqcdtools.base.logger as logger

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
    """ The following conversions take hc from PDG 2018. DOI: 10.1103/PhysRevD.98.030001. """
    if units=="MeVfm":
        return hcMeVfm
    elif units=="GeVfm":
        return hcGeVfm
    else:
        logger.TBError("Invalid unit specification.")
