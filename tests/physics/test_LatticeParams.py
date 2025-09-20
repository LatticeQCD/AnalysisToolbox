# 
# testLatticeParams.py                                                               
# 
# D. Clarke
# 
# Some tests of the lattice parameter class using 2021 scales.
#

from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.interfaces.HotQCD import HotQCDParams
from latqcdtools.interfaces.MILC import MILCParams
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def testLatticeParams():

    import params

    lpass = True

    lp = latticeParams(params.Ns, params.Nt, params.cbeta, params.cml, params.cms)
    lp.paramSummary()
    a = lp.geta()
    T = lp.getT()
    lpass *= print_results(a,0.14027137436021253,text='fK test, a')
    lpass *= print_results(T,175.84394864955703,text='fK test, T')
    del lp

    lp = HotQCDParams(params.Ns, params.Nt, params.cbeta, params.cml, params.cms, Nf='3')
    lp.paramSummary()
    logger.info('HotQCD:',lp.getcparams())
    del lp

    lp = MILCParams(params.Ns, params.Nt, 6.500, params.cml, params.cms, '411', Nf='211')
    logger.info('MILC:',lp.getcparams())
    del lp

    lp = MILCParams(24, 8, '35500','00239954','00499905','100481', Nf='111')
    logger.info('MILC Nf=1+1+1:',lp.getcparams())
    lp.paramSummary()
    del lp

    import params_r0

    lp = latticeParams(params_r0.Ns, params_r0.Nt, params_r0.cbeta, scaleType='r0', paramYear=2017)
    lp.paramSummary()
    a = lp.geta()
    T = lp.getT()
    lpass *= print_results(a,0.04195979097233082,text='r0 test, a')
    lpass *= print_results(T,146.96136335775122,text='r0 test, T')
    del lp

    concludeTest(lpass)


if __name__ == '__main__':
    testLatticeParams()