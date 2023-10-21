# 
# testJackknife.py                                                               
# 
# D. Clarke 
# 
# Quick test to make sure the jackknife works. Do not adjust the values of any of the variables, arrays, or arguments.
# 

from latqcdtools.statistics.statistics import std_err, std_mean
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.testing import print_results, concludeTest
import latqcdtools.base.logger as logger
import numpy as np

logger.set_log_level('INFO')

EPSILON=1e-15 # test precision

def simple_mean(a):
    return np.mean(a)

def div_simple(a, b, c):
    return b * a[0] / (c * a[1])

def div_old(a, b, c):
    return b * np.mean(a[0]) / (c * np.mean(a[1]))

def f(a, b, c):
    return b * np.mean(a[0]) / (c * np.mean(a[1]))

def div1(a, b, c):
    return ([[f(a, b, c), f(a, b, c)], [f(a, b, c), f(a, b, c)]],
            [f(a, b, c), f(a, b, c)], f(a, b, c))

def div2(a, b, c):
    return [[f(a, b, c), f(a, b, c)], [f(a, b, c), f(a, b, c)]]

def div3(a, b, c):
    return [f(a, b, c), f(a, b, c)]

def div4(a, b, c):
    return f(a, b, c)

def divnp1(a, b, c):
    return (np.array([[f(a, b, c), f(a, b, c)], [f(a, b, c), f(a, b, c)]]),
            np.array([f(a, b, c), f(a, b, c)]), np.array(f(a, b, c)))

A = np.array(range(1000))
B = np.array(range(1000, 2000))


testdata =  np.random.default_rng(234).normal(0, 1, 200)

def testJackknife():

    lpass = True

    # 1d jackknife test since bias is almost zero
    REFm = std_mean(testdata)
    REFe = std_err(testdata)
    TESTm, TESTe = jackknife(simple_mean, testdata, numb_blocks=1, conf_axis=0, nproc=1)
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "single proc 1 delete jackknife simple mean test", EPSILON) 

    REFm = 499.5
    REFe = 95.74271077563381
    TESTm, TESTe = jackknife(simple_mean, A, numb_blocks=10, conf_axis=0, nproc=1)
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "single proc simple mean test", EPSILON)
    TESTm, TESTe = jackknife(simple_mean, A, numb_blocks=10, conf_axis=0)
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "simple mean test", EPSILON)

    REFm = (np.array([[0.33583199323816093, 0.33583199323816093], [0.33583199323816093, 0.33583199323816093]]),
            np.array([0.33583199323816093, 0.33583199323816093]),
            0.33583199323816093323816093)
    REFe = (np.array([[0.042622475647127026, 0.042622475647127026], [0.042622475647127026, 0.042622475647127026]]),
            np.array([0.042622475647127026, 0.042622475647127026]),
            0.042622475647127026)

    TESTm, TESTe = jackknife(div1, [A, B], numb_blocks=10, args=(2, 2))
    lpass *= print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
                           TESTe[0].reshape(4), REFe[0].reshape(4), "div1, tuple[0]", EPSILON)
    lpass *= print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1, tuple[1]", EPSILON)

    lpass *= print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1, tuple[2]", EPSILON)
    TESTm, TESTe = jackknife(divnp1, [A, B], numb_blocks=10, args=(2, 2))
    lpass *= print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
                           TESTe[0].reshape(4), REFe[0].reshape(4), "div1np, tuple[0]", EPSILON)
    lpass *= print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1np, tuple[1]", EPSILON)
    lpass *= print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1np, tuple[2]", EPSILON)

    REFm = np.array([0.33583199323816093, 0.33583199323816093, 0.33583199323816093, 0.33583199323816093])
    REFe = np.array([0.042622475647127026, 0.042622475647127026, 0.042622475647127026, 0.042622475647127026])
    TESTm, TESTe = jackknife(div2, [A, B], numb_blocks=10, args=(2, 2))
    lpass *= print_results(TESTm.reshape(4), REFm, TESTe.reshape(4), REFe, "div2", EPSILON)

    REFm = np.array([0.33583199323816093, 0.33583199323816093])
    REFe = np.array([0.042622475647127026, 0.042622475647127026])
    TESTm, TESTe = jackknife(div3, [A, B], numb_blocks=10, args=(2, 2))
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "div3", EPSILON)

    REFm = 0.33583199323816093323816093
    REFe = 0.042622475647127026
    TESTm, TESTe = jackknife(div4, [A, B], numb_blocks=10, args=(2, 2))
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "div4", EPSILON)

    concludeTest(lpass)


if __name__ == '__main__':
    testJackknife()
