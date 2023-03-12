# 
# testUtilities.py                                                               
# 
# D. Clarke
# 
# Test some of the methods in the utilities module.
# 

import numpy as np
from latqcdtools.base.check import print_results
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import parallel_function_eval, naturalSort

logger.set_log_level('INFO')


def square(x):
    return x**2


def testUtilities():

    testArray = range(10)

    temp = parallel_function_eval(square,testArray,8)

    sum1 = 0.
    for i in testArray:
        sum1 += square(i)

    sum2 = np.sum(temp)

    print_results(sum1,sum2,text='parallel square test')

    testArray = ['thermalTable_mu0.0357', 'thermalTable_mu0.0952', 'thermalTable_mu0.1309', 'thermalTable_mu0.0833',
                 'thermalTable_mu0.0595', 'thermalTable_mu0.0119', 'thermalTable_mu0.0', 'thermalTable_mu0.0714',
                 'thermalTable_mu0.0476', 'thermalTable_mu0.1071', 'thermalTable_mu0.119', 'thermalTable_mu0.0238']
    sortArray = ['thermalTable_mu0.0', 'thermalTable_mu0.0119', 'thermalTable_mu0.119', 'thermalTable_mu0.0238',
                 'thermalTable_mu0.0357', 'thermalTable_mu0.0476', 'thermalTable_mu0.0595', 'thermalTable_mu0.0714',
                 'thermalTable_mu0.0833', 'thermalTable_mu0.0952', 'thermalTable_mu0.1071', 'thermalTable_mu0.1309']

    if naturalSort(testArray)==sortArray:
        logger.TBPass('natural sort')
    else:
        logger.TBFail('natural sort')


if __name__ == '__main__':
    testUtilities()