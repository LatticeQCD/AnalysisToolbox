# 
# testLegacy.py                                                               
# 
# Make sure legacy functionality doesn't break. 
# 
import numpy as np
from latqcdtools.testing import print_results, concludeTest
from latqcdtools.legacy import jackknife, bootstr, DEFAULTSEED


EPSILON=1e-15 # test precision


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


def testLegacyJack():

    lpass = True

    REFm = (np.array([[0.33583199323816093, 0.33583199323816093], [0.33583199323816093, 0.33583199323816093]]),
            np.array([0.33583199323816093, 0.33583199323816093]),
            0.33583199323816093323816093)
    REFe = (np.array([[0.042622475647127026, 0.042622475647127026], [0.042622475647127026, 0.042622475647127026]]),
            np.array([0.042622475647127026, 0.042622475647127026]),
            0.042622475647127026)
    REFsamp = [ (np.array([[0.13943564, 0.13943564],
                       [0.13943564, 0.13943564]]), np.array([0.13943564, 0.13943564]), 0.13943563633828227),
                (np.array([[0.18138663, 0.18138663],
                       [0.18138663, 0.18138663]]), np.array([0.18138663, 0.18138663]), 0.1813866331691294),
                (np.array([[0.22394803, 0.22394803],
                       [0.22394803, 0.22394803]]), np.array([0.22394803, 0.22394803]), 0.22394802608350695),
                (np.array([[0.26713323, 0.26713323],
                       [0.26713323, 0.26713323]]), np.array([0.26713323, 0.26713323]), 0.2671332348040387),
                (np.array([[0.31095608, 0.31095608],
                       [0.31095608, 0.31095608]]), np.array([0.31095608, 0.31095608]), 0.31095607533912917),
                (np.array([[0.35543077, 0.35543077],
                       [0.35543077, 0.35543077]]), np.array([0.35543077, 0.35543077]), 0.35543077471970763),
                (np.array([[0.40057199, 0.40057199],
                       [0.40057199, 0.40057199]]), np.array([0.40057199, 0.40057199]), 0.4005719863985231),
                (np.array([[0.44639481, 0.44639481],
                       [0.44639481, 0.44639481]]), np.array([0.44639481, 0.44639481]), 0.44639480634699735),
                (np.array([[0.49291479, 0.49291479],
                       [0.49291479, 0.49291479]]), np.array([0.49291479, 0.49291479]), 0.4929147898867918),
                (np.array([[0.54014797, 0.54014797],
                       [0.54014797, 0.54014797]]), np.array([0.54014797, 0.54014797]), 0.5401479692955027) ]


    TESTm, TESTe = jackknife(div1, [A, B], numb_blocks=10, args=(2, 2))
    lpass *= print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
                           TESTe[0].reshape(4), REFe[0].reshape(4), "div1, tuple[0]", EPSILON)
    lpass *= print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1, tuple[1]", EPSILON)

    lpass *= print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1, tuple[2]", EPSILON)
    samp, TESTm, TESTe = jackknife(divnp1, [A, B], numb_blocks=10, return_sample=True, args=(2, 2))
    lpass *= print_results(TESTm[0].reshape(4), REFm[0].reshape(4),
                           TESTe[0].reshape(4), REFe[0].reshape(4), "div1np, tuple[0]", EPSILON)
    lpass *= print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1np, tuple[1]", EPSILON)
    lpass *= print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1np, tuple[2]", EPSILON)

    for i in range(10):
        for j in range(len(samp[i])):
            lpass *= print_results(samp[i][j],REFsamp[i][j],text=f"div1np {i} {j}",prec=1e-7)

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


def testLegacyBoot():

    lpass = True

    REFm = ( np.array([[0.33364288472202835, 0.33364288472202835],[0.33364288472202835, 0.33364288472202835]]),
         np.array([0.33364288472202835, 0.33364288472202835]), 
         0.33364288472202847 )

    REFe = ( np.array([[0.006479415039739239, 0.006479415039739239],[0.006479415039739239, 0.006479415039739239]]),
         np.array([0.006479415039739239, 0.006479415039739239]), 
         0.006479415039739239 )

    TESTm, TESTe = bootstr(div1, [A, B], numb_samples=100, seed=DEFAULTSEED, args=(2, 2))
    lpass *= print_results(TESTm[0].reshape(4), REFm[0].reshape(4),TESTe[0].reshape(4), REFe[0].reshape(4), "div1, tuple[0]", EPSILON)
    lpass *= print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1, tuple[1]", EPSILON)
    lpass *= print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1, tuple[2]", EPSILON)
    TESTm, TESTe = bootstr(divnp1, [A, B], numb_samples=100, seed=DEFAULTSEED, args=(2, 2))
    lpass *= print_results(TESTm[0].reshape(4), REFm[0].reshape(4), TESTe[0].reshape(4), REFe[0].reshape(4), "div1np, tuple[0]", EPSILON)
    lpass *= print_results(TESTm[1], REFm[1], TESTe[1], REFe[1], "div1np, tuple[1]", EPSILON)
    lpass *= print_results(TESTm[2], REFm[2], TESTe[2], REFe[2], "div1np, tuple[2]", EPSILON)

    REFm = np.array([ 0.33364288472202835,  0.33364288472202835,  0.33364288472202835,  0.33364288472202835])
    REFe = np.array([0.006479415039739239, 0.006479415039739239, 0.006479415039739239, 0.006479415039739239])
    TESTm, TESTe = bootstr(div2, [A, B], numb_samples=100, seed=DEFAULTSEED, args=(2, 2))
    lpass *= print_results(TESTm.reshape(4), REFm, TESTe.reshape(4), REFe, "div2", EPSILON)

    REFm = np.array( [ 0.33364288472202835,  0.33364288472202835] )
    REFe = np.array( [0.006479415039739239, 0.006479415039739239] )
    TESTm, TESTe = bootstr(div3, [A, B], numb_samples=100, seed=DEFAULTSEED, args=(2, 2))
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "div3", EPSILON)

    REFm = 0.33364288472202847
    REFe = 0.006479415039739239
    TESTm, TESTe = bootstr(div4, [A, B], numb_samples=100, seed=DEFAULTSEED, args=(2, 2))
    lpass *= print_results(TESTm, REFm, TESTe, REFe, "div4", EPSILON)

    concludeTest(lpass)


if __name__ == '__main__':
    testLegacyJack()
    testLegacyBoot()