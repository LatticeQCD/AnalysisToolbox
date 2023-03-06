#
# testautocor.py
#
# D. Clarke
#
# Simple test for the getTauInt method to calculate integrated autocorrelation times. Compares against software
# by Bernd Berg.
#

from latqcdtools.statistics.statistics import getTauInt
from latqcdtools.base.check import print_results
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')


def testAutocor():

    nt = 48
    nbins = 8
    ts = [63.5, 62.1, 62.5, 63.5, 62.1, 62.1, 62.5, 62.5, 62.1, 62.5, 62.5, 62.5, 63.5, 63  , 63  , 63  , 63  , 63.5, 63  ,
          63.5, 63.5, 63.5, 63.5, 63.5, 62.5, 62.5, 62.5, 61.2, 62.1, 63.5, 63  , 63  , 62.5, 62.5, 62.5, 62.5, 63.5, 62.1,
          63.5, 63  , 62.1, 63  , 63.1, 63.7, 63.9, 62.4, 61.9, 61.6, 62  , 61.6, 62.9, 62.2, 62.7, 62.7, 62.2, 62.6, 62.1,
          62.1, 61.9, 62.6, 62.6, 62.4, 63.1, 62.5, 61.9, 62.2, 62.8, 61.8, 61.8, 62.1, 62.1, 61.5, 61.5, 62.1, 62.1, 62.3,
          61.9, 61.9, 61.9, 61.6, 61.9, 61.5, 61.9, 62.2, 61.8, 61.8, 61.8, 62.9, 62.1, 62.1, 62  , 61.9, 62.6, 62  , 62  ,
          62.4, 62.4, 62.2, 61.8, 62.3, 62  , 61.7, 61.7, 62  , 61.4, 61.4, 61.4, 62.3, 62  , 62  , 62  , 62.1]

    tau_int, tau_inte, tau_intbias, itpick = getTauInt(ts, nbins, nt, 'acor.d')

    TESTitpick   = 33
    TESTtau_int  = 18.24028851979112
    TESTtau_inte = 5.9003840732043
    TESTbias     = 4.79395572142051

    print_results( [TESTitpick,TESTtau_int,TESTtau_inte,TESTbias],[itpick,tau_int,tau_inte,tau_intbias], None, None, "tau_int", 1e-14 )


if __name__ == '__main__':
      testAutocor()