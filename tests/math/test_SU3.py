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
from latqcdtools.math.SU3 import SU3, id_3, ze_3
import latqcdtools.base.logger as logger
from latqcdtools.math.math import rel_check, invert 
from latqcdtools.testing import concludeTest


logger.set_log_level('INFO')

def testSU3():

    g = SU3(id_3)
    h = SU3(id_3)
    x = SU3(id_3)
    logger.debug('g=\n',g)
    logger.debug('h=\n',h)
    logger.debug('x=\n',x)

    lpass = True

    # Trace test
    if not rel_check(g.trace(),3):
        lpass = False
        logger.TBFail('Simple trace')

    # Assignment and comparison tests
    g[0,0]=2
    h[0,0]=2
    logger.debug('g=\n',g)

    if not rel_check(g,h): 
        lpass = False
        logger.TBFail('Assignment and comparison.')

    # Test mean along with some SU3 operations
    x[0,0]=(2+2+2*2+0+4)/5
    if not rel_check(  np.mean([g,h,2*g,g-h,h**2],axis=0) , x):
        lpass = False
        logger.TBFail('SU3 mean.')

    # Trace
    logger.debug('g=\n',g)
    if g.trace() != complex(4):
        lpass = False
        logger.TBFail('Trace:',g.trace())

    y = SU3([[ 0.17036869+0.89100163j,  0.03794601-0.13535497j,  0.28004516-0.098896j  ],
             [ 0.92906192-0.33496505j,  0.01763769+0.42305744j,  0.35911944+0.71267198j],
             [ 0.84188871-0.46038066j, -0.54672098+0.61514299j, -0.8191237 +0.77131243j] ])

    # Unitarize
    y.su3unitarize()
    if not y.isSU3():
        lpass = False
        logger.TBFail('su3unitarize')

    # Dagger
    if not rel_check(y.dagger(),invert(y)):
        lpass = False
        logger.TBFail('dagger')

    # Matrix power
    z=y**78
    if not z.isSU3():
        lpass=False
        logger.TBFail('matrix power')

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

    # Set to identity
    g.setToZero()
    if not rel_check(g,ze_3):
        lpass = False
        logger.TBFail('Set to zero.')

    concludeTest(lpass)


if __name__ == '__main__':
    testSU3()
