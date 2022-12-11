# 
# testWeightedMean.py                                                               
# 
# D. Clarke
# 
# Test of all the weighted mean methods.
#

import numpy as np
from numpy.random import normal
from latqcdtools.statistics.statistics import weighted_mean, weighted_mean_variance, unbiased_mean_variance, unbiased_sample_variance
from latqcdtools.base.check import print_results

PREC=1e-15

np.random.seed(7271978)
data = normal(10, 10, 100)
errors = np.full_like(data, 10)

weights = 1/errors**2

x         = weighted_mean(data, weights)
xe        = np.sqrt(weighted_mean_variance(errors))
xe_unbias = np.sqrt(unbiased_mean_variance(data, weights))
se_unbias = np.sqrt(unbiased_sample_variance(data, weights))

TRUEx        , TRUExe                  = 10.786789615738144, 1.0
TRUExe_unbias, TRUEse_unbias           = 1.1149172148742512, 11.14917214874251

print_results(TRUEx        ,x        ,text='weighted average',prec=PREC)
print_results(TRUExe       ,xe       ,text='weighted average error',prec=PREC)
print_results(TRUExe_unbias,xe_unbias,text='unbiased weighted average error',prec=PREC)
print_results(TRUEse_unbias,se_unbias,text='unbiased sample average error',prec=PREC)
