# 
# testUtilities.py                                                               
# 
# D. Clarke
# 
# Test some of the methods in the utilities module.
# 

import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import naturalSort, getMaxThreads 

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

    if naturalSort(testArray)==sortArray:
        logger.TBPass('natural sort')
    else:
        logger.TBFail('natural sort')

    logger.info('Testing getMaxThreads() functionality...')
    logger.info('Number of threads on system:',getMaxThreads())
    logger.TBPass('All tests passed')


if __name__ == '__main__':
    testUtilities()