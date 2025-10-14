# 
# testDeriv.py
# 
# H. Sandmeyer
# 
# Tests of methods to calculate derivatives numerically.
#

import numpy as np
from latqcdtools.math.num_deriv import diff_fit_grad, diff_fit_hess, diff_deriv
from latqcdtools.math.polynomials import Polynomial
from latqcdtools.testing import print_results, print_results_iter, concludeTest


def f(x, b):
    return b * x[0] * x[1] * np.exp(x[1])


def g(x, b):
    return [b * x[0] * x[1] * np.exp(x[1]), x[1]**2]


def grad_f(x, b):
    return [b * x[1] * np.exp(x[1]), b * x[0] * np.exp(x[1]) + b * x[0] * x[1] * np.exp(x[1])]


def jac_g(x, b):
    return [[b * x[1] * np.exp(x[1]), b * x[0] * np.exp(x[1]) + b * x[0] * x[1] * np.exp(x[1])], [0, 2 * x[1]]]


def hess_f(x, b):
    return b * np.array([[0, x[1] * np.exp(x[1]) + np.exp(x[1])], 
        [x[1] * np.exp(x[1]) + np.exp(x[1]), x[0] * (x[1] * np.exp(x[1]) + 2 * np.exp(x[1]))]])


def fit_func(x, params, b):
    return params[0]**2*x**2 + params[0]*params[1]*x + b/params[2]


def fit_grad_ref(x, params, b):
    return np.array([2*params[0]*x**2  + params[1]*x, params[0]*x, -b/params[2]**2])


def fit_hess_ref(x, params, b):
    return np.array([[2*x**2, x, 0], [x, 0, 0], [0, 0, 2*b/params[2]**3]])


def testDeriv():

    lpass = True

    test_param = [1, -2, 3]
    x   = 1
    opt = -1
    res = diff_fit_grad(x, test_param, fit_func, args=(opt,))
    res_true = fit_grad_ref(x, test_param, opt)
    lpass *= print_results_iter(res, res_true, "numerical gradient for fitting functions")

    res      = diff_fit_hess(x, test_param, fit_func, args=(opt,), h = 1e-5)
    res_true = fit_hess_ref(x, test_param, opt)
    lpass *= print_results_iter(res, res_true, "numerical hessian for fitting functions")

    x     = np.linspace(0,1,11)
    poly  = Polynomial([1 ,-1,1 ,-1, 1])
    dpoly = Polynomial([-1,2 ,-3, 4, 0])
    lpass *= print_results(dpoly(x), diff_deriv(x,poly), text="diff_deriv", prec=1e-6)

    concludeTest(lpass)


if __name__ == '__main__':
    testDeriv()
