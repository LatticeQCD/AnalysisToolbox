# 
# test_SU2.py                                                               
# 
# D. Clarke
# 

from latqcdtools.math.SU2 import SU2, sigma
from latqcdtools.math.math import rel_check, ze, exp
from latqcdtools.testing import concludeTest


def testSU2():

    g = SU2()
    g.setToRandom()

    lpass = True 

    theta = {}
    for i in [1,2,3]:
        theta[i] = g.projectPauli(i)
    
    arg = ze(2)
    for i in [1,2,3]:
        arg += theta[i]*sigma[i]
    
    recon = exp(1j*arg)

    lpass *= rel_check(g,recon)

    concludeTest(lpass)


if __name__ == '__main__':
    testSU2()
