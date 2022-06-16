# 
# test_deriv.py
# 
# H. Sandmeyer
# 
# Tests of methods to calculate derivatives numerically.
#


import numpy as np
from latqcdtools.math.num_deriv import diff_fit_grad, diff_fit_hess
from latqcdtools.base.check import rel_check
import latqcdtools.base.logger as logger

def print_results(res,  res_true, text):
    test = True
    res = np.asarray(res)
    res_true = np.asarray(res_true)
    it = np.nditer(res, flags=['multi_index'])
    while not it.finished:
        if not rel_check(res[it.multi_index], res_true[it.multi_index], 1e-4, 1e-6):
            test = False
            print("res[" + str(it.multi_index) + "] = " + str(res[it.multi_index])
                    + " != res_true[" + str(it.multi_index) + "] = " + str(res_true[it.multi_index]))
        it.iternext()
    if test:
        logger.TBPass(text)
    else:
        logger.TBFail(text)
    print()



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


print("\nTesting fitting gradient")
def fit_func(x, params, b):
    return params[0]**2*x**2 + params[0]*params[1]*x + b/params[2]

def func(params, b):
    return params[0]**2*x**2 + params[0]*params[1]*x + b/params[2]

def fit_grad_ref(x, params, b):
    return np.array([2*params[0]*x**2  + params[1]*x, params[0]*x, -b/params[2]**2])

def fit_hess_ref(x, params, b):
    return np.array([[2*x**2, x, 0], [x, 0, 0], [0, 0, 2*b/params[2]**3]])

test_param = [1, -2, 3]
x = 1
opt = -1
res = diff_fit_grad(x, test_param, fit_func, args=(opt,))
res_true = fit_grad_ref(x, test_param, opt)
print_results(res, res_true, "numerical gradient for fitting functions")

res = diff_fit_hess(x, test_param, fit_func, args=(opt,), eps = 1e-5)
res_true = fit_hess_ref(x, test_param, opt)
print_results(res, res_true, "numerical hessian for fitting functions")