#!/bin/python3

# 
# main_plotRatApprox.py                                                               
# 
# C. Schmidt 
# 
# A script for plotting rational approximation files made by SIMULATeQCD/src/tools/rational_approx/poly4.C.
# This gives a visual guide how the rational approximation compares to the exact result.
#

import sys
import numpy as np
from latqcdtools.base.plotting import plt, latexify, set_params
from latqcdtools.statistics.statistics import plot_func
import latqcdtools.base.logger as logger

if len(sys.argv) < 4:
    logger.TBError("Usage ", sys.argv[0], "in.rational m_l m_s")

ml = float(sys.argv[2])
ms = float(sys.argv[3])

Nmax = 15

r_inv_1f_num = np.zeros(Nmax)
r_inv_1f_den = np.zeros(Nmax)

r_1f_num     = np.zeros(Nmax)
r_1f_den     = np.zeros(Nmax)

r_bar_1f_num = np.zeros(Nmax)
r_bar_1f_den = np.zeros(Nmax)

r_inv_2f_num = np.zeros(Nmax)
r_inv_2f_den = np.zeros(Nmax)

r_2f_num     = np.zeros(Nmax)
r_2f_den     = np.zeros(Nmax)

r_bar_2f_num = np.zeros(Nmax)
r_bar_2f_den = np.zeros(Nmax)

# Read in rat_approx file, which holds the *_const variables
exec(open(sys.argv[1]).read())

def rat_func(x, r_const, num, den):
    ret = r_const
    for i in range(Nmax):
        ret += num[i]/(x+den[i])
    return ret

def diff_1f_rat_func(x, r_const, num, den, exp):
    return abs((x**exp - rat_func(x, r_const, num, den))/x**exp)

def exact_2f(x, exp):
    return ((x/(x+ms**2-ml**2))**exp)

def diff_2f_rat_func(x, r_const, num, den, exp):
    return abs((exact_2f(x, exp) - rat_func(x, r_const, num, den))/exact_2f(x,exp))

latexify()

set_params(ylogscale=True, xlabel="$x$", ylabel="(Approx - exact)/exact", title="$m_l = "+str(ml)+"$ $m_s = "+str(ms)+"$")

plot_func(diff_1f_rat_func, xmin = ms**2, xmax = 5, args=(r_inv_1f_const, r_inv_1f_num, r_inv_1f_den, 3/8.), label = "1f inv")
plot_func(diff_1f_rat_func, xmin = ms**2, xmax = 5, args=(r_1f_const, r_1f_num, r_1f_den, -3/8.), label = "1f")
plot_func(diff_1f_rat_func, xmin = ms**2, xmax = 5, args=(r_bar_1f_const, r_bar_1f_num, r_bar_1f_den, -3/4.), label = "1f bar")
plot_func(diff_2f_rat_func, xmin = ml**2, xmax = 5, args=(r_inv_2f_const, r_inv_2f_num, r_inv_2f_den, 1/4.), label = "2f inv")
plot_func(diff_2f_rat_func, xmin = ml**2, xmax = 5, args=(r_2f_const, r_2f_num, r_2f_den, -1/4.), label = "2f")
plot_func(diff_2f_rat_func, xmin = ms**2, xmax = 5, args=(r_bar_2f_const, r_bar_2f_num, r_bar_2f_den, -1/2.), label = "2f bar")

plt.show()
