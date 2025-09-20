# 
# MILC.py
# 
# D. Clarke
# 
# Some parameter combinations and naming conventions special to projects of the MILC collaboration.
#

from latqcdtools.interfaces.interfaces import HotQCD_MILC_Params

class MILCParams(HotQCD_MILC_Params):
    """
    A class to handle and check the input parameters of a lattice run, especially for MILC.
    """

    def __repr__(self) -> str:
        return "MILCParams"
