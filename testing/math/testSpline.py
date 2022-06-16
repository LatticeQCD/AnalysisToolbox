# 
# testSpline.py                                                               
# 
# D. Clarke
# 
# Test some of the convenience wrappers for spline methods.
# 

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LSQUnivariateSpline
from latqcdtools.math.spline import auto_knots,random_knots
from latqcdtools.base.plotting import plot_dots, plot_lines, set_params
from latqcdtools.base.check import print_results

x = np.linspace(-1, 1, 101)
y = 10*x**2 + np.random.randn(len(x))

knots = auto_knots(x, 3)

print_results(knots,[-0.5, 0.0, 0.5], text="auto_knots")

spline = LSQUnivariateSpline(x, y, knots, k=3)

plot_dots(x, y)
plot_lines(x,spline(x),marker=None)
set_params(xlabel='x',ylabel='y')

plt.savefig('autoSpline.pdf')
plt.clf()

knots=random_knots(x, 3, SEED=7271978)

print_results(knots,[-0.56, 0.11000000000000004, 0.5000000000000001], text="random_knots")

spline = LSQUnivariateSpline(x, y, knots, k=3)

plot_dots(x, y)
plot_lines(x,spline(x),marker=None)
set_params(xlabel='x',ylabel='y')
plt.savefig('randomSpline.pdf')
