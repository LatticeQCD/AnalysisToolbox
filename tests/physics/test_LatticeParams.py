# 
# test_LatticeParams.py                                                               
# 
# D. Clarke, K. Ebira
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

    # Issue #83 tests: get_aWorld
    from latqcdtools.physics.referenceScales import get_aWorld
    from latqcdtools.base.logger import ToolboxException
    from latqcdtools.physics.lattice_params import ignoreWorldWarning

    if get_aWorld('r0', 2012) != 'Nf21':
        logger.TBFail("get_aWorld('r0', 2012) != 'Nf21'")
        lpass = False
    if get_aWorld('r0', 2012.0) != 'Nf21':
        logger.TBFail("get_aWorld('r0', 2012.0) != 'Nf21'")
        lpass = False
    if get_aWorld('r0', 2017) != 'Nf0':
        logger.TBFail("get_aWorld('r0', 2017) != 'Nf0'")
        lpass = False
    if get_aWorld('r1', 2021) != 'Nf21':
        logger.TBFail("get_aWorld('r1', 2021) != 'Nf21'")
        lpass = False
    if get_aWorld('fk', 2021) != 'Nf21':
        logger.TBFail("get_aWorld('fk', 2021) != 'Nf21'")
        lpass = False
    if get_aWorld('t0') != 'Nf0':
        logger.TBFail("get_aWorld('t0') != 'Nf0'")
        lpass = False

    # Check error handling for unsupported scale or year
    try:
        get_aWorld('invalid_scale')
        logger.TBFail("get_aWorld failed to raise on invalid scaleType")
        lpass = False
    except ToolboxException:
        pass

    try:
        get_aWorld('r0', 1999)
        logger.TBFail("get_aWorld failed to raise on invalid year")
        lpass = False
    except ToolboxException:
        pass

    # Issue #83 tests: smart paramYear defaults and world disentanglement
    lp_r0_21 = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='r0', Nf='21')
    if lp_r0_21.year != 2012:
        logger.TBFail("Expected default paramYear 2012 for r0 with Nf='21', got", lp_r0_21.year)
        lpass = False
    if lp_r0_21.aWorld != 'Nf21':
        logger.TBFail("Expected aWorld 'Nf21', got", lp_r0_21.aWorld)
        lpass = False
    if lp_r0_21.physWorld != 'Nf21':
        logger.TBFail("Expected physWorld 'Nf21', got", lp_r0_21.physWorld)
        lpass = False
    if lp_r0_21.Nf != '21':
        logger.TBFail("Expected Nf '21', got", lp_r0_21.Nf)
        lpass = False
    del lp_r0_21

    # Issue #83 tests: pure gauge smart default
    lp_r0_quenched = latticeParams(Ns, Nt, 6.0, scaleType='r0', Nf=None)
    if lp_r0_quenched.year != 2017:
        logger.TBFail("Expected default paramYear 2017 for r0 with Nf=None, got", lp_r0_quenched.year)
        lpass = False
    if lp_r0_quenched.aWorld != 'Nf0':
        logger.TBFail("Expected aWorld 'Nf0', got", lp_r0_quenched.aWorld)
        lpass = False
    if lp_r0_quenched.physWorld != 'Nf21':
        logger.TBFail("Expected physWorld 'Nf21', got", lp_r0_quenched.physWorld)
        lpass = False
    if lp_r0_quenched.Nf is not None:
        logger.TBFail("Expected Nf None, got", lp_r0_quenched.Nf)
        lpass = False
    del lp_r0_quenched

    # Issue #83 tests: explicit world mismatch triggers warning
    lp_mismatch = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='r0', paramYear=2017, Nf='21')
    if lp_mismatch.aWorld != 'Nf0' or lp_mismatch.Nf != '21':
        logger.TBFail("Mismatch setup incorrect")
        lpass = False
    lp_mismatch.paramSummary()
    del lp_mismatch

    # Test scaleType='t0'
    lp_t0 = latticeParams(8, 4, 6.0, scaleType='t0')
    if lp_t0.physWorld != 'Nf21':
        logger.TBFail("Expected physWorld 'Nf21' for scaleType='t0', got", lp_t0.physWorld)
        lpass = False
    lp_t0.paramSummary()
    del lp_t0

    # Test ignoreWorldWarning()
    ignoreWorldWarning()
    lp_squelched = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='r0', paramYear=2017, Nf='21')
    lp_squelched.paramSummary()
    # Test zero-temperature ensemble (Nt=None) paramSummary execution
    lp_zero_temp = latticeParams(Ns, None, 6.5, scaleType='fk')
    lp_zero_temp.paramSummary()
    del lp_zero_temp

    concludeTest(lpass)


if __name__ == '__main__':
    testLatticeParams()
