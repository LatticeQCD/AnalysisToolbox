#!/usr/bin/env python3
from latqcdtools import spline_interpolate
from latqcdtools import plotting
from matplotlib import pyplot as plt

order = 3

xdata = [[5. , 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9],
        [5. , 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9],
        [5. , 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8]]

ydata = [[5. , 5.1, 5.1, 5.3, 5.4, 5.7, 5.6, 5.7, 5.8, 6.1],
        [6. , 6.1, 6.2, 6.3, 6.4, 6.6, 6.6, 6.7, 6.1, 6.9],
        [7. , 7.1, 7.1, 7.3, 7.0, 7.7, 7.8, 7.7, 7.7]]


def linear(x, a, b):
    return a*x + b


constraints = [[5.2, 1, 1], [5.2, 0, 6]]

knots, res, res_err, chi_dof = spline_interpolate.multi_spline_fit(linear, 2, [5, 6, 7], xdata, ydata,
        order = order, ncoeffs = None, knots = [5.3, 5.5, 5.6],
        tol = 1e-8, nsteady_deriv = 2, algorithms = ['curve_fit'], always_return = True)


plotting.plot_dots(xdata[0], ydata[0])
plotting.plot_dots(xdata[1], ydata[1])
plotting.plot_dots(xdata[2], ydata[2])

coeffs = [linear(5, *i) for i in res]
plotting.plot_func(spline_interpolate.constraint_spline, args = (knots, coeffs, order, 2, 0, constraints))


coeffs = [linear(6, *i) for i in res]
plotting.plot_func(spline_interpolate.constraint_spline, args = (knots, coeffs, order, 2, 0, constraints))


coeffs = [linear(7, *i) for i in res]
plotting.plot_func(spline_interpolate.constraint_spline, args = (knots, coeffs, order, 2, 0, constraints))

plt.show()
