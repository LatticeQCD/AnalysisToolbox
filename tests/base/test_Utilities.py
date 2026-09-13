# 
# test_Utilities.py                                                               
# 
# D. Clarke, K. Ebira
# 
# Test some of the methods in the utilities module.
# 

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.testing import concludeTest
from latqcdtools.base.logger import ToolboxException
from latqcdtools.base.utilities import comesBefore, naturalSort, envector, unvector, isArrayLike, \
    toNumpy, isFloatType, isComplexType, isScalar, isHigherDimensional, isIntType, isReal, \
    cleanOutput, printClean


def testUtilities():

    testArray = ['thermalTable_mu0.0357', 'thermalTable_mu0.0952', 'thermalTable_mu0.1309', 'thermalTable_mu0.0833',
                 'thermalTable_mu0.0595', 'thermalTable_mu0.0119', 'thermalTable_mu0.0'   , 'thermalTable_mu0.0714',
                 'thermalTable_mu0.0476', 'thermalTable_mu0.1071', 'thermalTable_mu0.119' , 'thermalTable_mu0.0238']
    sortArray = ['thermalTable_mu0.0'   , 'thermalTable_mu0.0119', 'thermalTable_mu0.0238', 'thermalTable_mu0.0357', 
                 'thermalTable_mu0.0476', 'thermalTable_mu0.0595', 'thermalTable_mu0.0714', 'thermalTable_mu0.0833', 
                 'thermalTable_mu0.0952', 'thermalTable_mu0.1071', 'thermalTable_mu0.119' , 'thermalTable_mu0.1309']

    lpass=True

    if naturalSort(testArray)!=sortArray:
        logger.TBFail('natural sort')
        for f in naturalSort(testArray):
            logger.TBFail(' ',f)
        lpass=False

    x = 3
    if x != unvector(envector(x)):
        logger.TBFail('unvector/envector')
        lpass=False

    if not isArrayLike(envector(x)):
        logger.TBFail('isArrayLike')
        lpass=False

    if not isReal(x):
        logger.TBFail('isReal')
        lpass=False

    if not isIntType(x):
        logger.TBFail('isIntType 1')
        lpass=False

    x = 3.143342342
    test = np.array([np.array([x])])
    if x != unvector(unvector(test)):
        logger.TBFail('unvector**2')
        lpass=False

    if isIntType(x):
        logger.TBFail('isIntType 2')
        lpass=False

    date1 = "2017/12/14 14:50:30"
    date2 = "2018/1/1 15:20:25"
    if not comesBefore(date1,date2):
        logger.TBFail('date comparison')
        lpass=False

    x1 = [1,1,1,1,1,1,1]
    x2 = [1,1,1,1,1]
    x3 = None 

    lpass *= not isHigherDimensional(x1)

    r1, r2, r3, r4 = toNumpy(x1,x2,x3,x2)

    if not type(r1) == np.ndarray:
        logger.TBFail('toNumpy r1')
        lpass=False
    if not type(r2) == np.ndarray:
        logger.TBFail('toNumpy r2')
        lpass=False
    if not r3 is None:
        logger.TBFail('toNumpy r3')
        lpass=False
    if not type(r4) == np.ndarray:
        logger.TBFail('toNumpy r4')
        lpass=False

    x1 = np.longdouble(1.)
    x2 = np.complex128(1.)
    lpass *= isFloatType(x1)
    lpass *= not isFloatType(x2)
    lpass *= isComplexType(x2)
    lpass *= not isComplexType(x1)
    lpass *= isScalar(x1)
    lpass *= isScalar(x2)
    lpass *= not isArrayLike(x1)
    lpass *= isReal(x1)
    lpass *= isReal(x2)

    # Test cleanOutput single integer sspace
    out1 = cleanOutput('hello', 'world', sspace=10)
    expected1 = '%10s  %10s' % ('hello', 'world')
    if out1 != expected1:
        logger.TBFail('cleanOutput sspace int: got', repr(out1), 'expected', repr(expected1))
        lpass = False

    # Test cleanOutput iterable sspace
    out2 = cleanOutput('alpha', 'beta', sspace=[8, 12])
    expected2 = '%8s  %12s' % ('alpha', 'beta')
    if out2 != expected2:
        logger.TBFail('cleanOutput sspace list: got', repr(out2), 'expected', repr(expected2))
        lpass = False

    # Test label spacing: label followed by column must not collide
    out_lbl = cleanOutput('col', label='LBL', sspace=5)
    expected_lbl = 'LBL  %5s' % 'col'
    if out_lbl != expected_lbl:
        logger.TBFail('cleanOutput label spacing: got', repr(out_lbl), 'expected', repr(expected_lbl))
        lpass = False

    # Test cleanOutput raises on array-like and non-scalar args
    for nested in ([1, 2], np.array([1, 2]), (1, 2), [], (), np.array([]),
                   range(0), range(3), set(), {1, 2}, {}, {'a': 1}, (x for x in range(3))):
        try:
            cleanOutput(nested)
            logger.TBFail('cleanOutput failed to raise on nested:', nested)
            lpass = False
        except ToolboxException:
            pass

    # Test isArrayLike and isHigherDimensional do not crash on dict/set
    if isArrayLike({}):
        logger.TBFail('isArrayLike({}) should be False')
        lpass = False
    if isArrayLike({'a': 1}):
        logger.TBFail('isArrayLike({"a": 1}) should be False')
        lpass = False
    if isArrayLike(set()):
        logger.TBFail('isArrayLike(set()) should be False')
        lpass = False
    if isHigherDimensional({}):
        logger.TBFail('isHigherDimensional({}) should be False')
        lpass = False

    # Test cleanOutput raises on insufficient iterable sspace
    try:
        cleanOutput('s1', 's2', sspace=[10])
        logger.TBFail('cleanOutput failed to raise on insufficient sspace')
        lpass = False
    except ToolboxException:
        pass

    # Test cleanOutput rejects bool as sspace
    try:
        cleanOutput('s1', sspace=True)
        logger.TBFail('cleanOutput failed to raise on sspace=True')
        lpass = False
    except ToolboxException:
        pass

    # Test isIntType and cleanOutput on unsigned numpy integers
    lpass *= isIntType(np.uint32(42))
    lpass *= not isIntType(True)
    lpass *= isScalar(np.uint32(42))
    out_uint = cleanOutput(np.uint32(42))
    if '4.20000000e+01' not in out_uint:
        logger.TBFail('cleanOutput failed to format np.uint32:', repr(out_uint))
        lpass = False

    # Test printClean executes cleanly
    try:
        printClean('test', 1.23, label='TEST_PRINT')
    except Exception as e:
        logger.TBFail('printClean raised unexpected exception:', e)
        lpass = False

    concludeTest(lpass)


if __name__ == '__main__':
    testUtilities()
