# 
# main_bootstrap.py                                                               
# 
# D. Clarke 
# 
# Quick example how flexible the bootstrap is. In this exercise, we have a set of 
# data and corresponding error bars. We are going to fit the data with splines, 
# using the Gaussian bootstrap to propagate error to some intermediate point.
#

import numpy as np
from latqcdtools.base.readWrite import readTable
from latqcdtools.math.spline import getSpline
from latqcdtools.base.plotting import plt, plot_dots, plot_fill, latexify
from latqcdtools.statistics.bootstr import bootstr_from_gauss
from latqcdtools.base.initialize import initialize, finalize

initialize()

latexify()

xdata, ydata, edata = readTable("../tests/statistics/wurf.dat", usecols=(0,2,3))

plot_dots(xdata,ydata,edata)

# We will interpolate to these x-values
xspline = np.linspace(np.min(xdata),np.max(xdata),101)

def splineError(data):
    """ We assume no errors in the xdata. For the ydata, we will pass the bootstrap
    routine ydata along with the errors. The input to this function, data, will then
    be generated in each bootstrap bin by drawing normally from ydata with a spread
    of edata. In this example we simply get a spline, but you wrap anything you want
    inside your bootstrap procedure.
    """
    ys = getSpline(xdata,data,3)
    return ys(xspline)

ybs, ybserr = bootstr_from_gauss(splineError,data=ydata,data_std_dev=edata,numb_samples=100)

plot_fill(xspline,ybs,ybserr,label='bootstrapped interpolation')

plt.show()

finalize()
