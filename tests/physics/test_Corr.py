# 
# test_Corr.py                                                               
# 
# D. Clarke
# 

from latqcdtools.physics.correlators import foldCorrelator
from latqcdtools.testing import concludeTest,print_results
from latqcdtools.base.cleanData import clipRange
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

# This is a possible way to handle when the correlator you get represents
# one where only every other timeslice was measured.
def simulate_dt(Nt):
    dtCorr = np.arange(0,Nt,2)
    toFold = []
    for i in range(Nt):
        if i%2==0:
            toFold.append(dtCorr[i//2])
        else:
            toFold.append(np.inf)
    toFold=np.array(toFold)
    return clipRange( foldCorrelator(toFold), maxVal=np.inf, allowEqual=False )


def testCorr():

    lpass = True

    test     = np.array([5,4,3,2,1,2,3,4])
    expected = np.array([5,4,3,2,1])
    rand     = rng.random(1000)

    lpass *= print_results(foldCorrelator(test),expected          ,text='basic test')
    lpass *= print_results(foldCorrelator(rand),foldSimple(rand)  ,text='test against simple implementation')
    lpass *= print_results(simulate_dt(10)     ,[0,5,5]           ,text='dt=2, Nt=10')
    lpass *= print_results(simulate_dt(20)     ,[0,10,10,10,10,10],text='dt=2, Nt=20')

    concludeTest(lpass)

if __name__=='__main__':
    testCorr()
