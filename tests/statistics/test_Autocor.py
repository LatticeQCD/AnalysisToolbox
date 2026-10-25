#
# testautocor.py
#
# D. Clarke
#
# Simple test for the getTauInt method to calculate integrated autocorrelation times. Compares against software
# by Bernd Berg.
#

import numpy as np
from latqcdtools.statistics.autocorrelation import getTauInt, gewekeZ
from latqcdtools.interfaces.simulationManagement import analyzeChain
from latqcdtools.base.initialize import TBRNG, DEFAULTSEED
import latqcdtools.base.logger as logger
from latqcdtools.testing import print_results, concludeTest


nt = 48
nbins = 8
ts = [63.5, 62.1, 62.5, 63.5, 62.1, 62.1, 62.5, 62.5, 62.1, 62.5, 62.5, 62.5, 63.5, 63  , 63  , 63  , 63  , 63.5, 63  ,
      63.5, 63.5, 63.5, 63.5, 63.5, 62.5, 62.5, 62.5, 61.2, 62.1, 63.5, 63  , 63  , 62.5, 62.5, 62.5, 62.5, 63.5, 62.1,
      63.5, 63  , 62.1, 63  , 63.1, 63.7, 63.9, 62.4, 61.9, 61.6, 62  , 61.6, 62.9, 62.2, 62.7, 62.7, 62.2, 62.6, 62.1,
      62.1, 61.9, 62.6, 62.6, 62.4, 63.1, 62.5, 61.9, 62.2, 62.8, 61.8, 61.8, 62.1, 62.1, 61.5, 61.5, 62.1, 62.1, 62.3,
      61.9, 61.9, 61.9, 61.6, 61.9, 61.5, 61.9, 62.2, 61.8, 61.8, 61.8, 62.9, 62.1, 62.1, 62  , 61.9, 62.6, 62  , 62  ,
      62.4, 62.4, 62.2, 61.8, 62.3, 62  , 61.7, 61.7, 62  , 61.4, 61.4, 61.4, 62.3, 62  , 62  , 62  , 62.1]
ts = np.array(ts)

def testAutocor():

    lpass = True

    tau_int, tau_inte, itpick = getTauInt(ts, nbins, nt, 'acor.d', showPlot=False)

    TESTitpick   = 33
    TESTtau_int  = 18.24028851979112
    TESTtau_inte = 5.7219342257755175

    lpass *= print_results( [TESTitpick,TESTtau_int,TESTtau_inte],[itpick,tau_int,tau_inte], None, None, "tau_int", 1e-14 )

    # Geweke test: an autocorrelated but stationary chain gives |z| of order 1, while a chain that starts far
    # from equilibrium and relaxes slowly gives a large |z|.
    rng      = TBRNG(DEFAULTSEED)
    nchain   = 5000
    chain    = np.zeros(nchain)
    for t in range(1,nchain):
        chain[t] = 0.9*chain[t-1] + rng.normal()
    drifting = chain + 5*np.exp(-np.arange(nchain)/1000)
    zstat    = gewekeZ(chain)
    zdrift   = gewekeZ(drifting)
    if abs(zstat) < 3 and abs(zdrift) > 3:
        logger.TBPass('Geweke test, z =',round(zstat,2),'for a stationary chain and',round(zdrift,2),'for a drifting one.')
    else:
        logger.TBFail('Geweke test, z =',zstat,'for a stationary chain and',zdrift,'for a drifting one.')
        lpass = False
    res = analyzeChain(np.arange(nchain), drifting, tpickMax=200)
    lpass *= print_results(res['gewekeZ'], zdrift, text='analyzeChain gewekeZ')

    concludeTest(lpass)


if __name__ == '__main__':
    testAutocor()
