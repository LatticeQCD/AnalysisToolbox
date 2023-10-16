# 
# testWeightedMean.py                                                               
# 
# D. Clarke
# 
# Test of all the weighted mean methods.
#

import numpy as np
from numpy.random import normal
from latqcdtools.statistics.statistics import gaudif, std_mean, std_err, weighted_mean, weighted_variance, \
    unbiased_mean_variance, unbiased_sample_variance
from latqcdtools.math.math import print_results
from latqcdtools.base.initialize import DEFAULTSEED
import latqcdtools.base.logger as logger


logger.set_log_level('INFO')
PREC=1e-15
np.random.seed(DEFAULTSEED)


def testWeightedMean():

    data   = normal(10, 10, 100)
    errors = np.full_like(data, 10)

    x         = weighted_mean(data, errors)
    xe        = np.sqrt(weighted_variance(errors))
    xe_unbias = np.sqrt(unbiased_mean_variance(data, errors))
    se_unbias = np.sqrt(unbiased_sample_variance(data, errors))

    TRUEx        , TRUExe        = 10.786789615738144, 1.0
    TRUExe_unbias, TRUEse_unbias = 1.1149172148742512, 11.14917214874251

    print_results(TRUEx        ,x        ,text='weighted average',prec=PREC)
    print_results(TRUExe       ,xe       ,text='weighted average error',prec=PREC)
    print_results(TRUExe_unbias,xe_unbias,text='unbiased weighted average error',prec=PREC)
    print_results(TRUEse_unbias,se_unbias,text='unbiased sample average error',prec=PREC)

    # Now see what happens we we combine two samples with different variances and statistics
    data1 = normal(0,1,600)
    data2 = normal(0,3,150)
    x1    = std_mean(data1)
    e1    = std_err(data1)
    x2    = std_mean(data2)
    e2    = std_err(data2)
    data  = np.concatenate((data1,data2))
    x     = std_mean(data)
    e     = std_err(data)
    x12   = weighted_mean([x1,x2],[e1,e2])
    e12   = np.sqrt(weighted_variance([e1,e2]))

    q = gaudif(x,e,x12,e12)

    if q<0.05 or e12>e:
        logger.TBFail('Combination test 1')
    else:
        logger.TBPass('Combination test 1') 

    N1  = len(data1)
    N2  = len(data2)
    w1  = 1/e1**2 / (1/e1**2 + 1/e2**2)
    w2  = 1/e2**2 / (1/e1**2 + 1/e2**2)
    w1 *= (N1+N2)/N1
    w2 *= (N1+N2)/N2

    data  = np.concatenate((data1*w1,data2*w2))
    x     = std_mean(data)
    e     = std_err(data)

    print_results(x,x12,e,e12,prec=1e-4,text='Combination test 2')

if __name__ == '__main__':
    testWeightedMean()