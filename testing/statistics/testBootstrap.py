# 
# testBootstrap.py
# 
# D. Clarke 
# 
# Quick test to make sure the bootstrap works. Do not adjust the values of any of the variables, arrays, or arguments.
# Do not adjust the SEED.
# 
from latqcdtools.statistics.bootstr import bootstr, bootstr_from_gauss
from latqcdtools.math.math import print_results
import numpy as np

EPSILON = 1e-16 # test precision

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

def div(a):
    return a[0]/a[1]

A, B = ( np.array(range(1000)),
         np.array(range(1000,2000)) )


def Test_Bootstrap():

    SEED = 196


    REFm = 498.69909
    REFe = 9.085239972364768
    TESTm, TESTe = bootstr(np.mean, A, numb_samples=100, seed=SEED, nproc=1)
    print_results(TESTm, REFm, TESTe, REFe, "single proc simple mean test", EPSILON)
    TESTm, TESTe = bootstr(np.mean, A, 100, seed=SEED)
    print_results(TESTm, REFm, TESTe, REFe, "simple mean test", EPSILON)

    REFm = ( np.array([[0.33272899963394864, 0.33272899963394864],[0.33272899963394864, 0.33272899963394864]]),
         np.array([0.33272899963394864, 0.33272899963394864]), 
         0.33272899963394875 )

    REFe = ( np.array([[0.00593633241664651, 0.00593633241664651],[0.00593633241664651, 0.00593633241664651]]),
         np.array([0.00593633241664651, 0.00593633241664651]), 
         0.00593633241664651 )

    TESTm, TESTe = bootstr(div1, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
    print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
              TESTe[0].reshape(4), REFe[0].reshape(4), "div1, tuple[0]", EPSILON)
    print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1, tuple[1]", EPSILON)
    print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1, tuple[2]", EPSILON)
    TESTm, TESTe = bootstr(divnp1, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
    print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
                   TESTe[0].reshape(4), REFe[0].reshape(4), "div1np, tuple[0]", EPSILON)
    print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1np, tuple[1]", EPSILON)
    print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1np, tuple[2]", EPSILON)

    REFm = np.array([0.33272899963394864, 0.33272899963394864, 0.33272899963394864, 0.33272899963394864])
    REFe = np.array([0.00593633241664651, 0.00593633241664651, 0.00593633241664651, 0.00593633241664651])
    TESTm, TESTe = bootstr(div2, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
    print_results(TESTm.reshape(4), REFm, TESTe.reshape(4), REFe, "div2", EPSILON)

    REFm = np.array( [0.33272899963394864, 0.33272899963394864] )
    REFe = np.array( [0.00593633241664651, 0.00593633241664651] )
    TESTm, TESTe = bootstr(div3, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
    print_results(TESTm, REFm, TESTe, REFe, "div3", EPSILON)

    REFm = 0.3327289996339486441208935
    REFe = 0.00593633241664651
    TESTm, TESTe = bootstr(div4, [A, B], numb_samples=100, seed=SEED, args=(2, 2))
    print_results(TESTm, REFm, TESTe, REFe, "div4", EPSILON)

    # Gaussian bootstrap tests

    TESTm, TESTe = bootstr_from_gauss(np.mean, data=[10], data_std_dev=[0.5], numb_samples=1000, seed=SEED)
    REFm = 9.992029809808589
    REFe = 0.4925667778426321
    print_results(TESTm, REFm, TESTe, REFe, "simple gauss", EPSILON)

    TESTm, TESTe = bootstr_from_gauss(div, data=[10,2], data_std_dev=[0.5,0.1], numb_samples=1000, seed=SEED)
    REFm = 5.005130586993748
    REFe = 0.3555674022912204
    print_results(TESTm, REFm, TESTe, REFe, "div gauss", EPSILON)


if __name__ == '__main__':
    Test_Bootstrap()