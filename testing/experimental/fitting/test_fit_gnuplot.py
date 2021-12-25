#!/usr/bin/env python3
from latqcdtools.fitting import *
from latqcdtools.readin import *
from latqcdtools.tools import print_results

print("\nComparing Toolbox fit with gnuplot:\n")

def fit_func(x, a, b, c):
    return a * x**2 + b * x + c

data = read_in("wurf.dat", 1, 3, 4)

res_gnu = fit_gnuplot(data, "a*x**2+b*x+c", ['a', 'b', 'c'], column_string="1:2:3 yerror",
        fit_range=[0,10])
res_fit = do_fit(fit_func, data[0], data[1], data[2])
print_results(res_gnu[0], res_fit[0], res_gnu[1], res_fit[1],
        "Comparing Gnuplot to own fit test")

res_gnu = fit_gnuplot(data[0:2], "a*x**2+b*x+c", ['a', 'b', 'c'], column_string="1:2",
        fit_range=[0,10])
res_fit = do_fit(fit_func, data[0], data[1])
print_results(res_gnu[0], res_fit[0], res_gnu[1], res_fit[1],
        "Comparing Gnuplot to own fit without errors")
