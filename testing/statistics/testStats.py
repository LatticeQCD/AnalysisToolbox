# 
# testStats.py                                                               
# 
# D. Clarke
# 
# Tests of some of the basic statistics methods.
#

import numpy as np
from latqcdtools.statistics.statistics import std_cov
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')


def testStats():

    A = []
    B = []
    for i in range(3):
        A.append(i)
        B.append(2-i)
    A = np.array(A)
    B = np.array(B)

    AB = np.vstack((A,B))

    C = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
           C[i][j] = np.mean( ( AB[i] - np.mean(AB[i]) )*( AB[j] - np.mean(AB[j]) ) )

    ltest = True

    compare = C==std_cov(AB)
    if not compare.all():
        logger.TBFail('std_cov')
        ltest = False

    if ltest:
        logger.TBPass('All tests passed.')


if __name__ == '__main__':
    testStats()