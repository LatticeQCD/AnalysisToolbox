# 
# test_LatticeParams.py                                                               
# 
# D. Clarke, K. Ebira, Google Antigravity
# 
# Some tests of the lattice parameter class using 2021 scales.
#

import os
import sys
import subprocess
import latqcdtools.physics.lattice_params as lp_mod
from latqcdtools.physics.lattice_params import latticeParams, ignoreWorldWarning
from latqcdtools.physics.referenceScales import get_aWorld, r0_div_a
from latqcdtools.physics.constants import fk_phys, r0_phys, convert
from latqcdtools.interfaces.collaborations import HotQCDParams, MILCParams
from latqcdtools.base.logger import ToolboxException
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

    # Issue #83 tests: get_aWorld and world handling

    if get_aWorld('r0', 2012) != 'Nf21':
        logger.TBFail("get_aWorld('r0', 2012) != 'Nf21'")
        lpass = False
    if get_aWorld('r0', '2012') != 'Nf21':
        logger.TBFail("get_aWorld('r0', '2012') != 'Nf21'")
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

    # Check error handling for unsupported scale or invalid year types (floats, NaN, inf, unknown years)
    try:
        get_aWorld('invalid_scale')
        logger.TBFail("get_aWorld failed to raise on invalid scaleType")
        lpass = False
    except ToolboxException:
        pass

    for bad_year in (1999, 2012.0, 2012.5, float('nan'), float('inf'), True):
        try:
            get_aWorld('r0', bad_year)
            logger.TBFail("get_aWorld failed to raise on invalid year:", bad_year)
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
    if lp_r0_quenched.world != 'SU3':
        logger.TBFail("Expected legacy world 'SU3' for pure gauge, got", lp_r0_quenched.world)
        lpass = False
    if lp_r0_quenched.physWorld != 'Nf21':
        logger.TBFail("Expected physWorld 'Nf21', got", lp_r0_quenched.physWorld)
        lpass = False
    if lp_r0_quenched.Nf is not None:
        logger.TBFail("Expected Nf None, got", lp_r0_quenched.Nf)
        lpass = False
    del lp_r0_quenched

    # Issue #83 tests: explicit world mismatch triggers warning, and ignoreWorldWarning() squelches it
    warn_messages = []
    orig_warn = logger.warn

    def _capture_warn(*args, **kwargs):
        warn_messages.append(" ".join(str(a) for a in args))
        return orig_warn(*args, **kwargs)

    orig_checkworld = lp_mod.CHECKWORLD
    try:
        logger.warn = _capture_warn
        lp_mod.CHECKWORLD = True
        lp_mismatch = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='r0', paramYear=2017, Nf='21')
        if lp_mismatch.aWorld != 'Nf0' or lp_mismatch.Nf != '21':
            logger.TBFail("Mismatch setup incorrect")
            lpass = False
        if not any("does not match parameterization world 'Nf0'" in m for m in warn_messages):
            logger.TBFail("Expected world-mismatch warning when CHECKWORLD=True")
            lpass = False
        lp_mismatch.paramSummary()
        del lp_mismatch

        # Invalid Nf must raise before emitting any world-mismatch warning
        warn_messages.clear()
        try:
            latticeParams(Ns, Nt, 6.5, scaleType='fk', Nf='invalid')
            logger.TBFail("Expected ToolboxException for invalid Nf")
            lpass = False
        except ToolboxException:
            pass
        if any("does not match parameterization world" in m for m in warn_messages):
            logger.TBFail("Unexpected world-mismatch warning emitted before Nf validation failed")
            lpass = False

        # Test ignoreWorldWarning() squelches mismatch warning without leaking global state
        ignoreWorldWarning()
        warn_messages.clear()
        lp_squelched = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='r0', paramYear=2017, Nf='21')
        if any("does not match parameterization world" in m for m in warn_messages):
            logger.TBFail("Expected ignoreWorldWarning() to suppress world-mismatch warning")
            lpass = False
        lp_squelched.paramSummary()
        del lp_squelched
    finally:
        logger.warn = orig_warn
        lp_mod.CHECKWORLD = orig_checkworld

    # Issue #83 tests: physWorld selector for fk ('nature', 'Nf211', 'Nf2', default 'Nf21')
    lp_fk_nature = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='fk', Nf='21', physWorld='nature')
    if lp_fk_nature.physWorld != 'nature':
        logger.TBFail("Expected physWorld 'nature', got", lp_fk_nature.physWorld)
        lpass = False
    lpass *= print_results(lp_fk_nature.fK, fk_phys(year=2018, units="MeV", world="nature"),
                           text='fK physWorld=nature default 2018')
    del lp_fk_nature

    lp_fk_2012 = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='fk', Nf='21', scaleYear=2012)
    if lp_fk_2012.physWorld != 'nature':
        logger.TBFail("Expected physWorld 'nature' for fk scaleYear=2012, got", lp_fk_2012.physWorld)
        lpass = False
    lpass *= print_results(lp_fk_2012.fK, fk_phys(year=2012, units="MeV", world="nature"),
                           text='fK scaleYear=2012 inferred nature')
    del lp_fk_2012

    lp_fk_2018 = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='fk', Nf='21', scaleYear=2018)
    if lp_fk_2018.physWorld != 'nature':
        logger.TBFail("Expected physWorld 'nature' for fk scaleYear=2018, got", lp_fk_2018.physWorld)
        lpass = False
    lpass *= print_results(lp_fk_2018.fK, fk_phys(year=2018, units="MeV", world="nature"),
                           text='fK scaleYear=2018 inferred nature')
    del lp_fk_2018

    # physWorld is not derived from Nf for fk when omitted, and can be selected explicitly
    lp_fk_211_def = latticeParams(Ns, Nt, 6.5, cml, cms, '411', scaleType='fk', Nf='211')
    if lp_fk_211_def.physWorld != 'Nf21':
        logger.TBFail("Expected default physWorld 'Nf21' for fk with Nf='211', got", lp_fk_211_def.physWorld)
        lpass = False
    del lp_fk_211_def

    lp_fk_211_exp = latticeParams(Ns, Nt, 6.5, cml, cms, '411', scaleType='fk', Nf='211', physWorld='Nf211')
    if lp_fk_211_exp.physWorld != 'Nf211':
        logger.TBFail("Expected explicit physWorld 'Nf211', got", lp_fk_211_exp.physWorld)
        lpass = False
    lpass *= print_results(lp_fk_211_exp.fK, fk_phys(year=2019, units="MeV", world="Nf211"),
                           text='fK explicit physWorld=Nf211')
    del lp_fk_211_exp

    # Degenerate Nf='2' works cleanly in latticeParams and with explicit physWorld='Nf2'
    lp_nf2 = latticeParams(Ns, Nt, 6.5, cml, cms, scaleType='fk', Nf='2', physWorld='Nf2')
    if lp_nf2.Nf != '2' or lp_nf2.physWorld != 'Nf2':
        logger.TBFail("Expected Nf='2' and physWorld='Nf2', got", lp_nf2.Nf, lp_nf2.physWorld)
        lpass = False
    lpass *= print_results(lp_nf2.fK, fk_phys(year=2019, units="MeV", world="Nf2"),
                           text='fK Nf=2 with physWorld=Nf2')
    del lp_nf2

    # r0 and r1 reject unsupported physWorld values
    for sc in ('r0', 'r1'):
        try:
            latticeParams(Ns, Nt, cbeta, cml, cms, scaleType=sc, Nf='21', physWorld='nature')
            logger.TBFail(f"Expected ToolboxException for scaleType='{sc}' with physWorld='nature'")
            lpass = False
        except ToolboxException:
            pass

    # Invalid scaleYear / paramYear types raise ToolboxException rather than raw ValueError/TypeError
    for bad_sy in ('bad', [2018], 2018.5, True):
        try:
            latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='fk', Nf='21', scaleYear=bad_sy)
            logger.TBFail("Expected ToolboxException for bad scaleYear:", bad_sy)
            lpass = False
        except ToolboxException:
            pass

    for bad_py in (2012.0, 2012.5, 'bad'):
        try:
            latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='r0', Nf='21', paramYear=bad_py)
            logger.TBFail("Expected ToolboxException for bad paramYear:", bad_py)
            lpass = False
        except ToolboxException:
            pass

    # Test scaleType='t0' keeps physWorld=None, rejects bogus scaleYear, raises cleanly in geta/getT, and keeps paramSummary safe
    try:
        latticeParams(8, 4, 6.0, scaleType='t0', scaleYear=9999)
        logger.TBFail("Expected ToolboxException for scaleType='t0' with invalid scaleYear=9999")
        lpass = False
    except ToolboxException:
        pass

    lp_t0 = latticeParams(8, 4, 6.0, scaleType='t0')
    if lp_t0.physWorld is not None:
        logger.TBFail("Expected physWorld None for scaleType='t0', got", lp_t0.physWorld)
        lpass = False
    try:
        lp_t0.geta()
        logger.TBFail("Expected ToolboxException from lp_t0.geta()")
        lpass = False
    except ToolboxException:
        pass
    try:
        lp_t0.getT()
        logger.TBFail("Expected ToolboxException from lp_t0.getT()")
        lpass = False
    except ToolboxException:
        pass
    lp_t0.paramSummary()
    del lp_t0

    # Numerical scale checks for r0 (Nf='21' vs pure gauge) and r1
    r0_fm = r0_phys(year=2014, units='fm')
    lp_r0_21_num = latticeParams(24, 8, 6.285, scaleType='r0', Nf='21')
    a_r0_21_exp = r0_fm / r0_div_a(6.285, 2012)
    lpass *= print_results(lp_r0_21_num.geta(), a_r0_21_exp, text='r0 Nf=21 (2012) a')
    lpass *= print_results(lp_r0_21_num.getT(), convert(1.0 / (a_r0_21_exp * 8), 'fminv', 'MeV'), text='r0 Nf=21 (2012) T')
    del lp_r0_21_num

    lp_r0_q_num = latticeParams(24, 8, 6.285, scaleType='r0', Nf=None)
    a_r0_q_exp = r0_fm / r0_div_a(6.285, 2017)
    lpass *= print_results(lp_r0_q_num.geta(), a_r0_q_exp, text='r0 pure-gauge (2017) a')
    lpass *= print_results(lp_r0_q_num.getT(), convert(1.0 / (a_r0_q_exp * 8), 'fminv', 'MeV'), text='r0 pure-gauge (2017) T')
    del lp_r0_q_num

    lp_r1 = latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='r1', Nf='21')
    if lp_r1.aWorld != 'Nf21' or lp_r1.physWorld != 'Nf21' or lp_r1.scaleYear != 2010:
        logger.TBFail("Unexpected r1 worlds/scaleYear:", lp_r1.aWorld, lp_r1.physWorld, lp_r1.scaleYear)
        lpass = False
    if lp_r1.geta() <= 0 or lp_r1.getT() <= 0:
        logger.TBFail("Invalid r1 a/T:", lp_r1.geta(), lp_r1.getT())
        lpass = False
    del lp_r1

    # Invalid physWorld on fk and conflicting physWorld/scaleYear pairs raise ToolboxException
    try:
        latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='fk', Nf='21', physWorld='bogus')
        logger.TBFail("Expected ToolboxException for invalid physWorld on fk")
        lpass = False
    except ToolboxException:
        pass

    try:
        latticeParams(Ns, Nt, cbeta, cml, cms, scaleType='fk', Nf='21', physWorld='nature', scaleYear=2019)
        logger.TBFail("Expected ToolboxException for conflicting physWorld='nature' and scaleYear=2019")
        lpass = False
    except ToolboxException:
        pass

    # Test zero-temperature ensemble (Nt=None) paramSummary execution
    lp_zero_temp = latticeParams(Ns, None, 6.5, scaleType='fk')
    lp_zero_temp.paramSummary()
    del lp_zero_temp

    # CLI smoke tests for applications/main_getBeta.py and applications/main_getTempAndSpacing.py
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    app_beta = os.path.join(repo_root, 'applications', 'main_getBeta.py')
    app_temp = os.path.join(repo_root, 'applications', 'main_getTempAndSpacing.py')
    cli_env = os.environ.copy()
    cli_env.pop('DISPLAY', None)

    res_beta_r0 = subprocess.run([sys.executable, app_beta, '--Nt', '8', '--scale', 'r0', '--T', '300'],
                                 capture_output=True, text=True, env=cli_env)
    if res_beta_r0.returncode != 0 or 'beta_sol = 6.038' not in res_beta_r0.stdout or 'aWorld = Nf0' not in res_beta_r0.stdout:
        logger.TBFail("main_getBeta.py --scale r0 default compat failed:", res_beta_r0.stdout, res_beta_r0.stderr)
        lpass = False

    res_beta_fk = subprocess.run([sys.executable, app_beta, '--Nt', '8', '--scale', 'fk', '--T', '155', '--physWorld', 'nature'],
                                 capture_output=True, text=True, env=cli_env)
    if res_beta_fk.returncode != 0 or 'physWorld = nature' not in res_beta_fk.stdout or 'scaleYear = 2018' not in res_beta_fk.stdout:
        logger.TBFail("main_getBeta.py --scale fk --physWorld nature failed:", res_beta_fk.stdout, res_beta_fk.stderr)
        lpass = False

    res_beta_bad = subprocess.run([sys.executable, app_beta, '--Nt', '8', '--scale', 'r0', '--T', '300', '--Nf', 'bogus'],
                                  capture_output=True, text=True, env=cli_env)
    if res_beta_bad.returncode == 0:
        logger.TBFail("main_getBeta.py failed to reject --Nf bogus")
        lpass = False

    res_temp_r0 = subprocess.run([sys.executable, app_temp, '--Nt', '8', '--beta', '6.285', '--scale', 'r0'],
                                 capture_output=True, text=True, env=cli_env)
    if res_temp_r0.returncode != 0 or '437.37' not in res_temp_r0.stdout or 'Nf0 (r0 2017)' not in res_temp_r0.stdout:
        logger.TBFail("main_getTempAndSpacing.py --scale r0 default compat failed:", res_temp_r0.stdout, res_temp_r0.stderr)
        lpass = False

    concludeTest(lpass)


if __name__ == '__main__':
    testLatticeParams()
