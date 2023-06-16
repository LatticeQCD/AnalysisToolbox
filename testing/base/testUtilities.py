# 
# testUtilities.py                                                               
# 
# D. Clarke
# 
# Test some of the methods in the utilities module.
# 

import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import naturalSort, getMaxThreads, envector, unvector, isArrayLike 

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

    ltest=True

    if naturalSort(testArray)!=sortArray:
        logger.TBFail('natural sort')
        ltest=False

    logger.info('Testing getMaxThreads() functionality...')
    logger.info('Number of threads on system:',getMaxThreads())

    x = 3
    if x != unvector(envector(x)):
        logger.TBFail('unvector/envector')
        ltest=False

    if not isArrayLike(envector(x)):
        logger.TBFail('isArrayLike')
        ltest=False

    if ltest:
        logger.TBPass('All tests passed.')
    else:
        logger.TBError('At least one test failed.')


if __name__ == '__main__':
    testUtilities()