# 
# testStats.py                                                               
# 
# D. Clarke
# 
# Tests of some of the basic statistics methods.
#

import numpy as np
from latqcdtools.statistics.statistics import gaudif, studif
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

eps=1e-7

# Some measurements and error
x1=0.4655
e1=0.0088
x2=0.501
e2=0.045
x3=0.480
e3=0.023
x4=0.4745
e4=0.0063

# Some invented numbers of data
ndat1 = 15
ndat2 = 9

# Results produced by software of "Markov Chain Monte Carlo Simulations and Their Statistical Analysis, World
# Scientific, Singapore, 2004.
q12control=0.4387984
q13control=0.5559897
q14control=0.4056413
s12control=0.33726853


def testStats():

    ltest = True

    A = []
    B = []
    for i in range(3):
        A.append(i)
        B.append(2-i)
    A = np.array(A)
    B = np.array(B)

    q12=gaudif(x1,e1,x2,e2)
    q13=gaudif(x1,e1,x3,e3)
    q14=gaudif(x1,e1,x4,e4)

    # Test gaussian difference
    if abs(q12-q12control) > eps:
        ltest=False
    if abs(q13-q13control) > eps:
        ltest=False
    if abs(q14-q14control) > eps:
        ltest=False

    s12=studif(x1,e1,ndat1,x2,e2,ndat2)

    # Test student difference
    if abs(s12-s12control) > eps:
        ltest=False

    # Student and Gaussian difference tests should agree for large sample sizes
    if studif(x1,e1,300,x2,e2,300)/gaudif(x1,e1,x2,e2) > 1.005:
        ltest=False
        

    if ltest:
        logger.TBPass('All tests passed.')
    else:
        logger.TBError('At least one test failed.')


if __name__ == '__main__':
    testStats()