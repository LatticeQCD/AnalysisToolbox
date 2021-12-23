#!/usr/bin/env python3
from numpy import *
from latqcdtools.readin import *
from latqcdtools.spline_interpolate import *
import ctypes
from latqcdtools.plotting import *
import matplotlib.pyplot as plt

xdata, ydata, ydata_err=np.array(read_in("test.dat"))
xdata=np.array(xdata)
ydata=np.array(ydata)
ydata_err=np.array(ydata_err)

xspline, res=bspline_sci(xdata, ydata, ydata_err, npoints = 1000)

plot_dots(xdata, ydata, ydata_err)
plot_lines(xspline, res, marker=None)
plt.show()

xspline, res, res_err=bspline(xdata, ydata, ydata_err, npoints = 1000)

plt.clf()
plot_dots(xdata, ydata, ydata_err)
plot_fill(xspline, res, res_err)
plt.show()

xspline, res, res_err=bspline(xdata, ydata, ydata_err, xspline = [1, 2, 2.25, 2.5, 3])

plt.clf()
plot_dots(xdata, ydata, ydata_err)
plot_fill(xspline, res, res_err)
plt.show()
