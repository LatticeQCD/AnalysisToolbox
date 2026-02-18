# 
# collaborations.py                                                               
# 
# D. Clarke
# 
# Some convenience methods when analyzing data from various lattice collaborations. They generically
# have some conventions, for example when naming ensembles. These methods let you quickly generate
# e.g. ensemble names following those conventions. 
# 

from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.base.utilities import substringBetween
from latqcdtools.base.check import checkType
import latqcdtools.base.logger as logger

class HotQCD_MILC_Params(latticeParams):
    """ 
    A class to handle and check the input parameters of a lattice run using conventions common to both the
    HotQCD and MILC collaborations. 
    """

    def __repr__(self) -> str:
        return "HotQCD_MILC_Params"

    # String often used to label lattice configurations.
    def getcgeom(self):
        return 'l'+str(self.Ns)+str(self.Nt)
    def getcparams(self):
        if self.Nf=='211' or self.Nf=='111':
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2+'m'+self.cm3
        else:
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2


def paramFromEnsLabel(ensemble,format='MILC'):
    """ 
    Given an ensemble string, get the parameters out of it. 

    Args:
        ensemble (str): ensemble label

    Returns:
        tuple: Ns, Nt, Nf, beta string, mass1 string, mass2 string, mass3 string
    """
    checkType(str,ensemble=ensemble)
    checkType(str,format=format)
    Ns, Nt, Nf, cbeta, cm1, cm2, cm3 = None, None, None, None, None, None, None
    if format=='MILC':
        NsNt = substringBetween(ensemble,'l','f') 
        if len(NsNt)==3:
            Ns=NsNt[:2]
            Nt=NsNt[-1]
        elif len(NsNt)==4:
            Ns=NsNt[:2]
            Nt=NsNt[2:]
        else:
            logger.TBRaise('I assume Ns has 2 digits and Nt has 1-2 digits.') 
        Nf    = substringBetween(ensemble,'f','b') 
        cbeta = substringBetween(ensemble,'b','m') 
        cm1   = ensemble.split('m')[1].strip()
        if len(Nf)>1:
            cm2 = ensemble.split('m')[2].strip()
        if len(Nf)>2:
            cm3 = ensemble.split('m')[3].strip()
    else:
        logger.TBRaise('Unsupported format',format)
    return int(Ns), int(Nt), Nf, cbeta, cm1, cm2, cm3 


class HotQCDParams(HotQCD_MILC_Params):
    """
    A class to handle and check the input parameters of a lattice run, especially for HotQCD.
    """
    def __repr__(self) -> str:
        return "HotQCDParams"
    def getcGradFlowPureGauge(self) -> str:
        return 's'+str(self.Ns).zfill(3)+'t'+str(self.Nt).zfill(2)+'_b0'+self.cbeta+'00'


class MILCParams(HotQCD_MILC_Params):
    """
    A class to handle and check the input parameters of a lattice run, especially for MILC.
    """

    def __repr__(self) -> str:
        return "MILCParams"
