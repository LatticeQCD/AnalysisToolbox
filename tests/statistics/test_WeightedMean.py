# 
# testWeightedMean.py                                                               
# 
# D. Clarke
# 
# Test of all the weighted mean methods.
#

import numpy as np
from latqcdtools.statistics.statistics import gaudif, std_mean, std_err, weighted_mean, weighted_variance, \
    unbiased_mean_variance, unbiased_sample_variance
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.base.initialize import DEFAULTSEED, TBRNG
import latqcdtools.base.logger as logger


PREC=1e-15


def testWeightedMean():

    rng = TBRNG(DEFAULTSEED)

    lpass = True

    data   = rng.normal(10, 10, 100)
    errors = np.full_like(data, 10)

    x         = weighted_mean(data, errors)
    xe        = np.sqrt(weighted_variance(errors))
    xe_unbias = np.sqrt(unbiased_mean_variance(data, errors))
    se_unbias = np.sqrt(unbiased_sample_variance(data, errors))

    TRUEx        , TRUExe        = 10.679990952078827, 1.0
    TRUExe_unbias, TRUEse_unbias = 0.9274155524028469, 9.274155524028467

    lpass *= print_results(TRUEx        ,x        ,text='weighted average',prec=PREC)
    lpass *= print_results(TRUExe       ,xe       ,text='weighted average error',prec=PREC)
    lpass *= print_results(TRUExe_unbias,xe_unbias,text='unbiased weighted average error',prec=PREC)
    lpass *= print_results(TRUEse_unbias,se_unbias,text='unbiased sample average error',prec=PREC)

    # Unequal errors, compared against the reliability-weights formula written out explicitly.
    data   = rng.normal(5, 2, 8)
    errors = np.array([0.5,1,1,2,3,4,6,10])
    w      = 1/errors**2
    V1     = np.sum(w)
    V2     = np.sum(w**2)
    mu     = w @ data / V1
    s2     = w @ (data-mu)**2 / (V1 - V2/V1)
    lpass *= print_results(s2        ,unbiased_sample_variance(data, errors),text='unbiased sample variance, unequal errors',prec=1e-14)
    lpass *= print_results(s2*V2/V1**2,unbiased_mean_variance(data, errors) ,text='unbiased mean variance, unequal errors'  ,prec=1e-14)

    # Now see what happens we we combine two samples with different variances and statistics
    data1 = rng.normal(0,1,600)
    data2 = rng.normal(0,3,150)
    x1    = std_mean(data1)
    e1    = std_err(data1)
    x2    = std_mean(data2)
    e2    = std_err(data2)
    data  = np.concatenate((data1,data2))
    x     = std_mean(data)
    e     = std_err(data)
    x12   = weighted_mean(np.array([x1,x2]),np.array([e1,e2]))
    e12   = np.sqrt(weighted_variance(np.array([e1,e2])))

    q = gaudif(x,e,x12,e12)

    if q<0.05 or e12>e:
        logger.TBFail('Combination test 1')
        lpass = False

    N1  = len(data1)
    N2  = len(data2)
    w1  = 1/e1**2 / (1/e1**2 + 1/e2**2)
    w2  = 1/e2**2 / (1/e1**2 + 1/e2**2)
    w1 *= (N1+N2)/N1
    w2 *= (N1+N2)/N2

    data  = np.concatenate((data1*w1,data2*w2))
    x     = std_mean(data)
    e     = std_err(data)

    lpass *= print_results(x,x12,e,e12,prec=1e-3,text='Combination test 2')

    concludeTest(lpass)


if __name__ == '__main__':
    testWeightedMean()
