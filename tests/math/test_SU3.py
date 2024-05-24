# 
# testSU3.py                                                               
# 
# D. Clarke
# 
# Test some of the basic methods for SU3 matrices.
# 

import numpy as np
from latqcdtools.base.speedify import numbaON
numbaON()
from latqcdtools.math.SU3 import SU3, id_3
import latqcdtools.base.logger as logger
from latqcdtools.math.math import rel_check
from latqcdtools.testing import concludeTest


logger.set_log_level('INFO')

def testSU3():

    g = SU3(id_3)
    h = SU3(id_3)
    x = SU3(id_3)

    lpass = True

    # Assignment and comparison tests
    g[0,0]=2
    h[0,0]=2

    if not rel_check(g,h): 
        lpass = False
        logger.TBFail('Assignment and comparison.')

    # Test mean along with some SU3 operations
    x[0,0]=(2+2+2*2+0+4)/5
    if not rel_check(  np.mean([g,h,2*g,g-h,h**2],axis=0) , x):
        lpass = False
        logger.TBFail('SU3 mean.')

    # Trace
    if g.trace() != complex(4):
        lpass = False
        logger.TBFail('Trace:',g.trace())

    # Making a random matrix
    g.setToRandom()
    if not g.isSU3():
        lpass = False
        logger.TBFail('Set to random.')

    # Set to identity
    g.setToIdentity()
    if not rel_check(g,id_3):
        lpass = False
        logger.TBFail('Set to identity.')

    concludeTest(lpass)


if __name__ == '__main__':
    testSU3()
