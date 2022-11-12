# 
# interfaces.py
# 
# D. Clarke
# 
# Some common classes and functions that may be shared among multiple interfaces modules.
#


import latqcdtools.base.logger as logger
from latqcdtools.physics.lattice_params import latticeParams


class HotQCD_MILC_Params(latticeParams):
    """ A class to handle and check the input parameters of a lattice run using conventions common to both the
        HotQCD and MILC collaborations. """

    # String often used to label lattice configurations.
    def getcgeom(self):
        return 'l'+str(self.Ns)+str(self.Nt)
    def getcparams(self):
        if self.Nf=='211':
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2+'m'+self.cm3
        else:
            return self.getcgeom()+'f'+str(self.Nf)+'b'+self.cbeta+'m'+self.cm1+'m'+self.cm2
