# 
# testJackknife.py                                                               
# 
# D. Clarke 
# 
# Quick test to make sure the jackknife works. Do not adjust the values of any of the variables, arrays, or arguments.
# 

from latqcdtools.statistics.statistics import std_err, std_mean, binSeries
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger
import numpy as np

logger.set_log_level('INFO')

EPSILON=1e-15 # test precision

A = np.array(range(1000))


def testJackknife():

    testdata = np.random.default_rng(234).normal(0, 1, 200)

    lpass = True

    # 1d jackknife test since bias is almost zero
    REFm = std_mean(testdata)
    REFe = std_err(testdata)
    TESTm, TESTe = jackknife(std_mean, testdata, numb_blocks=len(testdata), conf_axis=0) 
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "delete-1 jackknife simple mean test", EPSILON) 

    REFm = 499.5
    REFe = 95.74271077563381
    REFsamp = [49.5, 149.5, 249.5, 349.5, 449.5, 549.5, 649.5, 749.5, 849.5, 949.5]
    samp, TESTm, TESTe = jackknife(std_mean, A, numb_blocks=10, conf_axis=0, return_sample=True)
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "simple mean test", EPSILON)
    lpass *= print_results(samp,REFsamp,text='simple mean sample',prec=EPSILON)

    TESTm, TESTe = jackknife(std_mean, binSeries(A,10), numb_blocks=10, conf_axis=0) 
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "remove-N vs N-bin with remove-1", EPSILON)

    concludeTest(lpass)


if __name__ == '__main__':
    testJackknife()
