#!/usr/bin/env python3

#
# testErrorProp.py                                                               
# 
# H. Sandmeyer
# 
# Quick test of the error propagation function.
#

import matplotlib.pyplot as plt
import numpy as np
from latqcdtools.base.plotting import error_prop_func, plot_func
from latqcdtools.base.check import print_results
from latqcdtools.statistics.statistics import error_prop


def func(x, A, B, OPT):
    return OPT*A**2*np.sin(x)+B

def grad(x, A, B, OPT):
    return [OPT*2*A*np.sin(x), 1]

def err_func(x, A, B, A_err, B_err, OPT):
    return np.sqrt((OPT*2*A*np.sin(x)*A_err)**2+B_err**2)


x_test = [0.5, 0.1, 2]
a = 1
b = 2
a_err = 0.1
b_err = 0.2
opt = 2

res_true = err_func(x_test, 1, 2, a_err, b_err, opt)

res = error_prop_func(x_test, func, [1, 2], [a_err, b_err], args=(opt,))
print_results(res, res_true, text = "Error_prop using diff quotient")

res = error_prop(lambda p, OPT: func(x_test, *p, OPT), [1, 2], [a_err, b_err], args=(opt,))[1]
print_results(res, res_true, text = "Direct error propagation using lambda to wrap x")

print_results(res, res_true, text = "Error_prop using self made grad")
res_true = err_func(x_test, 1, 2, a_err, b_err, opt)

plot_func(func, args = [a, b, opt], func_err = err_func, args_err=[a,b,a_err,b_err, opt])
plot_func(func, args = [a, b, opt], args_err = [a_err,b_err])
plot_func(func, args = [a, b, opt], args_err = [a_err,b_err], func_sup_numpy = True)
plot_func(func, args = [a, b, opt], args_err = [a_err,b_err,opt], grad = grad,
        func_sup_numpy = True, title = "Please check if all error bands are the same")

plt.show()