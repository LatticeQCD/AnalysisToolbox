# 
# test_Corr.py                                                               
# 
# D. Clarke
# 

from latqcdtools.physics.correlators import foldCorrelator
from latqcdtools.testing import concludeTest,print_results
from latqcdtools.base.initialize import rng
import numpy as np

def foldSimple(corr):
    N = len(corr)
    res = []
    res.append(corr[0])
    for i in range(1,int(N/2)):
        re = corr[i]
        rep = corr[N-i]
        res.append((re+rep)/2.)
    res.append(corr[int(N/2)])
    return np.array(res)


def testCorr():

    lpass = True

    test     = np.array([5,4,3,2,1,2,3,4])
    expected = np.array([5,4,3,2,1])

    lpass *= print_results(foldCorrelator(test),expected,text='basic test')

    rand = rng.random(1000)

    lpass *= print_results(foldCorrelator(rand),foldSimple(rand),text='test against simple implementation')

    concludeTest(lpass)

if __name__=='__main__':
    testCorr()
