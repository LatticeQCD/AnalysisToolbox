# 
# MILC.py
# 
# D. Clarke
# 
# Some parameter combinations and naming conventions special to projects of the MILC code.
#

from latqcdtools.interfaces.interfaces import HotQCD_MILC_Params

class MILCParams(HotQCD_MILC_Params):
    """A class to handle and check the input parameters of a lattice run, especially for MILC."""

    def __init__(self, Nsigma, Ntau, coupling, mass1=None, mass2=None, mass3=None, scaleType='fk', paramYear=2021, Nf='21', scaleYear=2019):
        HotQCD_MILC_Params.__init__(self, Nsigma, Ntau, coupling, mass1, mass2, mass3, scaleType, paramYear, Nf, scaleYear)
        if not isinstance(coupling,str):
            self.cbeta = str(int(coupling*100))

    def __repr__(self) -> str:
        return "MILCParams"
