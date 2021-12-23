# 
# testJackknife.py                                                               
# 
# D. Clarke 
# 
# Quick test to make sure the jackknife works. Do not adjust the values of any of the
# variables, arrays, or arguments. 
# 
from latqcdtools.jackknife import jackknife, jackknife_old
from latqcdtools.tools import print_results
import numpy as np


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


a, b, c = ( np.array(range(1000)), 
            np.array(range(1000,2000)),
            np.array(range(2000,3000)) )


REFm = 499.5
REFe = 95.74271077563381
TESTm, TESTe = jackknife(simple_mean, a, numb_blocks=10, conf_axis=1, args=(), cov=False, parallelize=False)
print_results(TESTm, REFm, TESTe, REFe, "single proc simple mean test", 1e-16)
TESTm, TESTe = jackknife(simple_mean, a, numb_blocks=10, conf_axis=1) 
print_results(TESTm, REFm, TESTe, REFe, "simple mean test", 1e-16)


REFm = 0.33583199323816093
REFe = 0.042622475647127026
TESTm, TESTe = jackknife_old(div_old, [a, b], numb_blocks=10, args=(2, 2))
print_results(TESTm, REFm, TESTe, REFe, "jackknife_old test", 1e-16)


REFm = ( np.array([[0.33583199, 0.33583199],[0.33583199, 0.33583199]]), 
         np.array([0.33583199, 0.33583199]), 
         0.33583199323816093 )
REFe = ( np.array([[0.04262248, 0.04262248],[0.04262248, 0.04262248]]), 
         np.array([0.04262248, 0.04262248]), 
         0.042622475647127026 ) 
TESTm, TESTe = jackknife(div1, [a, b], numb_blocks=10, args=(2, 2))
print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
                   TESTe[0].reshape(4), REFe[0].reshape(4), "div1, tuple[0]", 1e-6) 
print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1, tuple[1]", 1e-6)
print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1, tuple[2]", 1e-6)
TESTm, TESTe = jackknife(divnp1, [a, b], numb_blocks=10, args=(2, 2))
print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
                   TESTe[0].reshape(4), REFe[0].reshape(4), "div1np, tuple[0]", 1e-6) 
print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1np, tuple[1]", 1e-6)
print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1np, tuple[2]", 1e-6)


REFm = np.array([0.33583199, 0.33583199, 0.33583199, 0.33583199])
REFe = np.array([0.04262248, 0.04262248, 0.04262248, 0.04262248])
TESTm, TESTe = jackknife(div2, [a, b], numb_blocks=10, args=(2, 2))
print_results(TESTm.reshape(4), REFm, TESTe.reshape(4), REFe, "div2", 1e-6)


REFm = np.array( [0.33583199, 0.33583199] )
REFe = np.array( [0.04262248, 0.04262248] )
TESTm, TESTe = jackknife(div3, [a, b], numb_blocks=10, args=(2, 2))
print_results(TESTm, REFm, TESTe, REFe, "div3", 1e-6)

REFm = 0.33583199323816093
REFe = 0.042622475647127026
TESTm, TESTe = jackknife(div4, [a, b], numb_blocks=10, args=(2, 2))
print_results(TESTm, REFm, TESTe, REFe, "div4", 1e-6)
