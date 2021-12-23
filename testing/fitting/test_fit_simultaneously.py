#!/usr/bin/env python3
from latqcdtools.simultaneous_fit import *
from latqcdtools.statistics import *
from latqcdtools.plotting import *
import numpy as np

xdata_1 = np.linspace(0, 2*np.pi, 100)
xdata_2 = np.linspace(np.pi, 3*np.pi, 100)

raw_data_1 = np.random.normal(loc = 2 * np.cos(1 * xdata_1),
        scale = 1, size = (1000, len(xdata_1)))
raw_data_2 = np.random.normal(loc = 2 * np.cos(2 * xdata_1),
        scale = 1, size = (1000, len(xdata_2)))

ydata_1, edata_1 = mean_and_err(raw_data_1)
ydata_2, edata_2 = mean_and_err(raw_data_2)

xdata = [xdata_1, xdata_2]
ydata = [ydata_1, ydata_2]
edata = [edata_1, edata_2]


def fit_func(x, add_param, params):
    if add_param == 1:
        return params[0] * np.cos(params[1] * x)
    else:
        return params[0] * np.cos(params[2] * x)


res, res_err, chi_dof = simultaneous_fit(fit_func,
        [1, 0], xdata, ydata, edata, expand = False, start_params = [2, 1, 2])

print("First example:")
print(res, res_err)

plot_dots(xdata_1, ydata_1, edata_1)
plot_dots(xdata_2, ydata_2, edata_2)

plot_lines(xdata_1, fit_func(xdata_1, 1, res), marker = 'None')
plot_lines(xdata_2, fit_func(xdata_2, 0, res), marker = 'None')

plt.show()




print("Second example:")
def linear(x, a, b):
    return x*a + b


def fit_func(x, add_param, params):
    a = linear(add_param, params[0], params[1])
    b = linear(add_param, params[2], params[3])
    return  a * np.cos(b * x)

res, res_err, chi_dof = simultaneous_fit(fit_func,
        [2, 1], xdata, ydata, edata, expand = False,
        start_params = [0, 2, -1, 3])


print(res, res_err)
plt.clf()

plot_dots(xdata_1, ydata_1, edata_1)
plot_dots(xdata_2, ydata_2, edata_2)

plot_lines(xdata_1, fit_func(xdata_1, 2, res), marker = 'None')
plot_lines(xdata_2, fit_func(xdata_2, 1, res), marker = 'None')

plt.show()

