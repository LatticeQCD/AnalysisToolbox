# 
# testBootstrap.py
# 
# D. Clarke 
# 
# Quick test to make sure the bootstrap works. Do not adjust the values of any of the variables, arrays, or arguments.
# Do not adjust the SEED.
# 
from latqcdtools.statistics.bootstr import bootstr
from latqcdtools.base.check import print_results
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


A, B = ( np.array(range(1000)),
         np.array(range(1000,2000)) )

SEED = 196 


REFm = 498.69909
REFe = 9.085239972364768
TESTm, TESTe = bootstr(np.mean, A, numb_samples=100, seed=SEED, parallelize=False)
print_results(TESTm, REFm, TESTe, REFe, "single proc simple mean test", 1e-16)
TESTm, TESTe = bootstr(np.mean, A, 100, seed=SEED)
print_results(TESTm, REFm, TESTe, REFe, "simple mean test", 1e-16)


REFm = ( np.array([[0.33272899963394864, 0.33272899963394864],[0.33272899963394864, 0.33272899963394864]]), 
         np.array([0.33272899963394864, 0.33272899963394864]), 
         0.33272899963394875 ) 
REFe = ( np.array([[0.00593633241664651, 0.00593633241664651],[0.00593633241664651, 0.00593633241664651]]), 
         np.array([0.00593633241664651, 0.00593633241664651]), 
         0.00593633241664651 )
TESTm, TESTe = bootstr(div1, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
              TESTe[0].reshape(4), REFe[0].reshape(4), "div1, tuple[0]", 1e-6) 
print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1, tuple[1]", 1e-6)
print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1, tuple[2]", 1e-6)
TESTm, TESTe = bootstr(divnp1, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
                   TESTe[0].reshape(4), REFe[0].reshape(4), "div1np, tuple[0]", 1e-6) 
print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1np, tuple[1]", 1e-6)
print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1np, tuple[2]", 1e-6)


REFm = np.array([0.33272899963394864, 0.33272899963394864, 0.33272899963394864, 0.33272899963394864])
REFe = np.array([0.00593633241664651, 0.00593633241664651, 0.00593633241664651, 0.00593633241664651])
TESTm, TESTe = bootstr(div2, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
print_results(TESTm.reshape(4), REFm, TESTe.reshape(4), REFe, "div2", 1e-6)


REFm = np.array( [0.33272899963394864, 0.33272899963394864] )
REFe = np.array( [0.00593633241664651, 0.00593633241664651] )
TESTm, TESTe = bootstr(div3, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
print_results(TESTm, REFm, TESTe, REFe, "div3", 1e-6)

REFm = 0.3327289996339486441208935 
REFe = 0.00593633241664651
TESTm, TESTe = bootstr(div4, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
print_results(TESTm, REFm, TESTe, REFe, "div4", 1e-6)