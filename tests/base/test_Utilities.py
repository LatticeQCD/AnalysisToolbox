# 
# testUtilities.py                                                               
# 
# D. Clarke
# 
# Test some of the methods in the utilities module.
# 

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.testing import concludeTest
from latqcdtools.base.utilities import comesBefore, naturalSort, envector, unvector, isArrayLike, \
    toNumpy, isFloatType, isComplexType, isScalar, isHigherDimensional

logger.set_log_level('INFO')


def square(x):
    return x**2


def testUtilities():

    testArray = ['thermalTable_mu0.0357', 'thermalTable_mu0.0952', 'thermalTable_mu0.1309', 'thermalTable_mu0.0833',
                 'thermalTable_mu0.0595', 'thermalTable_mu0.0119', 'thermalTable_mu0.0'   , 'thermalTable_mu0.0714',
                 'thermalTable_mu0.0476', 'thermalTable_mu0.1071', 'thermalTable_mu0.119' , 'thermalTable_mu0.0238']
    sortArray = ['thermalTable_mu0.0'   , 'thermalTable_mu0.0119', 'thermalTable_mu0.119' , 'thermalTable_mu0.0238',
                 'thermalTable_mu0.0357', 'thermalTable_mu0.0476', 'thermalTable_mu0.0595', 'thermalTable_mu0.0714',
                 'thermalTable_mu0.0833', 'thermalTable_mu0.0952', 'thermalTable_mu0.1071', 'thermalTable_mu0.1309']

    lpass=True

    if naturalSort(testArray)!=sortArray:
        logger.TBFail('natural sort')
        lpass=False

    x = 3
    if x != unvector(envector(x)):
        logger.TBFail('unvector/envector')
        lpass=False

    if not isArrayLike(envector(x)):
        logger.TBFail('isArrayLike')
        lpass=False

    x = 3.143342342
    test = np.array([np.array([x])])
    if x != unvector(unvector(test)):
        logger.TBFail('unvector**2')
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

    x1 = np.float128(1.)
    x2 = np.complex128(1.)
    lpass *= isFloatType(x1)
    lpass *= not isFloatType(x2)
    lpass *= isComplexType(x2)
    lpass *= not isComplexType(x1)
    lpass *= isScalar(x1)
    lpass *= isScalar(x2)
    lpass *= not isArrayLike(x1)

    concludeTest(lpass)


if __name__ == '__main__':
    testUtilities()