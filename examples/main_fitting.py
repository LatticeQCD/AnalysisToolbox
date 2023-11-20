# 
# fitExample.py
# 
# H. Sandmeyer
# 
# Example of a simple 3-parameter quadratic fit.
#


import numpy as np
import matplotlib.pyplot as plt
from latqcdtools.statistics.fitting import Fitter
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.initialize import initialize, finalize
import latqcdtools.base.logger as logger


initialize('example_fitting.log')


# Here we define our fit function. we pass it its independent variable followed by the fit parameters we are
# trying to determine.
def fit_func(x, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a*x**2 + b*x + c


xdata, ydata = readTable("../testing/statistics/wurf.dat", usecols=(0,2)) 

# We initialize our Fitter object. The function has to look like
#            func(x, params, *args).
fitter = Fitter(fit_func, xdata, ydata)

# Here we try a fit, using the 'curve_fit' method, specifying the starting guesses for the fit parameters. Since
# detailedInfo = True, we will get back additional statistical information like the log of the Gaussian Bayes factor. 
res, res_err, chi_dof, stats = fitter.try_fit(start_params = [1, 2, 3], algorithms = ['curve_fit'], detailedInfo = True)

logger.info(" a , b,  c : ",res)
logger.info(" ae, be, ce: ",res_err)
logger.info("chi2/d.o.f.: ",chi_dof)
logger.info("     logGBF: ",stats['logGBF'])
logger.info("       pcov: \n\n",stats['pcov'],"\n")

# We can plot the fit and data using commands like this one. You can combine these commands with anything from
# plotting.py if you want to spice up your plot a bit.
fitter.plot_fit(domain=(np.min(xdata),np.max(xdata)))
fitter.plot_data()
plt.show()

finalize()


