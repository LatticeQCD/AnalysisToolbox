# 
# testBootstrap.py
# 
# D. Clarke 
# 
# Quick test to make sure the bootstrap works. Do not adjust the values of any of the variables, arrays, or arguments.
# 

from latqcdtools.statistics.bootstr import bootstr, bootstr_from_gauss
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.base.initialize import DEFAULTSEED
import numpy as np

EPSILON = 1e-16 # test precision

def simple_mean(a):
    return np.mean(a)

def div(a):
    return a[0]/a[1]

A =  np.array(range(1000))


def Test_Bootstrap():

    lpass = True

    REFm =  500.3124500000001
    REFe =  9.202595107612115
    TESTm, TESTe = bootstr(np.mean, A, numb_samples=100, seed=DEFAULTSEED, nproc=1)
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "single proc simple mean test", EPSILON)
    TESTm, TESTe = bootstr(np.mean, A, 100, seed=DEFAULTSEED)
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "simple mean test", EPSILON)

    # Gaussian bootstrap tests
    TESTm, TESTe = bootstr_from_gauss(np.mean, data=[10], data_std_dev=[0.5], numb_samples=1000, err_by_dist=False, seed=DEFAULTSEED)
    REFm = 9.994645407149935
    REFe =  0.4946243316482945
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "simple gauss", EPSILON)

    TESTm, TESTe = bootstr_from_gauss(div, data=[10,2], data_std_dev=[0.5,0.1], numb_samples=1000, err_by_dist=False, seed=DEFAULTSEED)
    REFm = 5.004629053835201
    REFe = 0.3474936630487199
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "div gauss", EPSILON)

    concludeTest(lpass)


if __name__ == '__main__':
    Test_Bootstrap()