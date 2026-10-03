# 
# testBootstrap.py
# 
# D. Clarke 
# 
# Quick test to make sure the bootstrap works. Do not adjust the values of any of the variables, arrays, or arguments.
# 

from latqcdtools.statistics.bootstr import bootstr, bootstr_from_gauss
from latqcdtools.statistics.statistics import KSTest_1side
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.base.initialize import DEFAULTSEED
import latqcdtools.base.logger as logger
import numpy as np
import scipy as sp

EPSILON = 1e-16 # test precision

def simple_mean(a):
    return np.mean(a)

def div(a):
    return a[0]/a[1]

A =  np.array(range(1000))


def testBootstrap():

    lpass = True

    # Test that nothing changes
    REFm =  498.91077000000007
    REFe =  9.048125857459329
    samp, TESTm, TESTe = bootstr(np.mean, A, numb_samples=100, seed=DEFAULTSEED, nproc=1, return_sample=True, err_by_dist=False)
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "single proc simple mean test", EPSILON)

    # Test that the bootstrap distribution is reasonable
    normalCDF = sp.stats.norm(loc=REFm,scale=REFe).cdf
    if KSTest_1side(samp,normalCDF) < 0.05:
        lpass = False
        logger.TBFail('Significant KS tension')

    TESTm, TESTe = bootstr(np.mean, A, 100, seed=DEFAULTSEED, err_by_dist=False)
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "simple mean test", EPSILON)

    # Gaussian bootstrap tests
    TESTm, TESTe = bootstr_from_gauss(np.mean, data=[10], data_std_dev=[0.5], numb_samples=1000, err_by_dist=False, seed=DEFAULTSEED)
    REFm = 9.9943152721107
    REFe =  0.49608800937424563
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "simple gauss", EPSILON)

    TESTm, TESTe = bootstr_from_gauss(div, data=[10,2], data_std_dev=[0.5,0.1], numb_samples=1000, err_by_dist=False, seed=DEFAULTSEED)
    REFm = 4.978630190450123
    REFe = 0.35750848724496614
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "div gauss", EPSILON)

    # same_rand_for_obs=True should make the observables perfectly correlated
    samp, _, _ = bootstr_from_gauss(lambda x: x, data=[10,2], data_std_dev=[0.5,0.1], numb_samples=1000,
                                    same_rand_for_obs=True, seed=DEFAULTSEED, return_sample=True)
    samp = np.array(samp)
    if not np.isclose(np.corrcoef(samp[:,0],samp[:,1])[0,1],1):
        lpass = False
        logger.TBFail('same_rand_for_obs gauss not perfectly correlated')

    # Six observables with identical data, resampled along conf_axis=2
    B = np.tile(A, (2,3,1))
    for same_rand in [True, False]:
        samp, _, _ = bootstr(lambda x: np.mean(x,axis=-1).ravel(), B, 100, conf_axis=2, same_rand_for_obs=same_rand,
                             seed=DEFAULTSEED, return_sample=True)
        samp = np.array(samp)
        if np.allclose(samp, samp[:,[0]]) != same_rand:
            lpass = False
            logger.TBFail('same_rand_for_obs =',same_rand,'gives wrong correlation for conf_axis=2')

    concludeTest(lpass)


if __name__ == '__main__':
    testBootstrap()
