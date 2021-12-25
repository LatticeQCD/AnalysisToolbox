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


def fit_func(x, add_param, params, index_array, default_array):
    my_params = list(default_array)
    for i, index in enumerate(index_array):
        my_params[index] = params[i]
    if add_param == 1:
        return my_params[0] * np.cos(my_params[1] * x)
    else:
        return my_params[0] * np.cos(my_params[2] * x)


index_array = [2,]
default_array = [2, 1, 2]
res, res_err, chi_dof = simultaneous_fit(fit_func,
        [1, 0], xdata, ydata, edata, expand = False, start_params = [2,], args = (index_array, default_array))

print(res, res_err)

plot_dots(xdata_1, ydata_1, edata_1)
plot_dots(xdata_2, ydata_2, edata_2)

plot_lines(xdata_1, fit_func(xdata_1, 1, res, index_array, default_array), marker = 'None')
plot_lines(xdata_2, fit_func(xdata_2, 0, res, index_array, default_array), marker = 'None')

plt.show()



