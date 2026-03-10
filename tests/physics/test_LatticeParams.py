# 
# testLatticepy                                                               
# 
# D. Clarke
# 
# Some tests of the lattice parameter class using 2021 scales.
#

from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.interfaces.collaborations import HotQCDParams, MILCParams
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')

Ns=32
Nt=8
cbeta='6500'
cml='00347'
cms='09369'

def testLatticeParams():

    lpass = True

    lp = latticeParams(Ns, Nt, cbeta, cml, cms, Nf='21')
    lp.paramSummary()
    a = lp.geta()
    T = lp.getT()
    lpass *= print_results(a,0.14027137436021253,text='fK test, a')
    lpass *= print_results(T,175.84394864955703,text='fK test, T')
    del lp

    lp = HotQCDParams(Ns, Nt, cbeta, cml, cms, Nf='3')
    lp.paramSummary()
    logger.info('HotQCD:',lp.getcparams())
    del lp

    lp = MILCParams(Ns, Nt, 6.500, cml, cms, '411', Nf='211')
    logger.info('MILC:',lp.getcparams())
    del lp

    lp = MILCParams(24, 8, '35500','00239954','00499905','100481', Nf='111')
    logger.info('MILC Nf=1+1+1:',lp.getcparams())
    lp.paramSummary()
    del lp

    concludeTest(lpass)


if __name__ == '__main__':
    testLatticeParams()
