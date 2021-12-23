#!/usr/bin/env python3

from latqcdtools.tools import *
from latqcdtools.plotting import *
from latqcdtools.statistics import error_prop
import matplotlib.pyplot as plt
import numpy as np
from latqcdtools.num_deriv import algopy_avail


def func(x, a, b, opt):
    return opt*a**2*np.sin(x)+b

def grad(x, a, b, opt):
    return [opt*2*a*np.sin(x), 1]

def err_func(x, a, b, a_err, b_err, opt):
    return np.sqrt((opt*2*a*np.sin(x)*a_err)**2+b_err**2)


x_test = [0.5, 0.1, 2]
a = 1
b = 2
a_err = 0.1
b_err = 0.2
opt = 2

res_true = err_func(x_test, 1, 2, a_err, b_err, opt)

res = error_prop_func(x_test, func, [1, 2], [a_err, b_err], args=(opt,))
print_results(res, res_true, text = "Error_prop using diff quotient")

res = error_prop(lambda p, opt: func(x_test, *p, opt), [1, 2], [a_err, b_err],
    args=(opt,), use_diff = True)[1]
print_results(res, res_true, text = "Direct error propagation using lambda to wrap x")

print_results(res, res_true, text = "Error_prop using self made grad")
res_true = err_func(x_test, 1, 2, a_err, b_err, opt)


plot_func(func, args = [a, b, opt], func_err = err_func, args_err=[a,b,a_err,b_err, opt])
plot_func(func, args = [a, b, opt], args_err = [a_err,b_err])
plot_func(func, args = [a, b, opt], args_err = [a_err,b_err], func_sup_numpy = True)
plot_func(func, args = [a, b, opt], args_err = [a_err,b_err,opt], grad = grad,
        func_sup_numpy = True, title = "Please check if all error bands are the same")

if algopy_avail:
    res = error_prop_func(x_test, func, [1, 2], [a_err, b_err], args=(opt,), use_diff = False)
    print_results(res, res_true, text = "Error_prop using algopy")
    res = error_prop_func(x_test, func, [1, 2], [a_err, b_err], grad=grad, args=(opt,))
    plot_func(func, args=[a, b, opt], func_err=error_prop_func,
            args_err=[func, [a,b], [a_err,b_err], None, False, (opt,)])

plt.show()