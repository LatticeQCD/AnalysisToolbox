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
import latqcdtools.base.logger as logger

logger.set_log_level('INFO')

print("\n Example of a simple 3-parameter quadratic fit.\n")

# Here we define our fit function. we pass it its independent variable followed by the fit parameters we are
# trying to determine.
def fit_func(x, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a*x**2 + b*x + c


def fitExample():

    xdata, ydata, edata = np.genfromtxt("wurf.dat", usecols=(0,2,3), unpack=True)

    # We initialize our Fitter object. If expand = True, fit_func has to look like
    #            func(x, a, b, *args)
    #        otherwise it has to look like
    #            func(x, params, *args).
    fitter = Fitter(fit_func, xdata, ydata, expand = False)

    # Here we try a fit, using the 'curve_fit' method, specifying the starting guesses for the fit parameters. Since
    # ret_pcov = True, we will get back the covariance matrix as well.
    res, res_err, chi_dof, pcov = fitter.try_fit(start_params = [1, 2, 3], algorithms = ['curve_fit'], ret_pcov = True)

    print(" a , b,  c : ",res)
    print(" ae, be, ce: ",res_err)
    print("chi2/d.o.f.: ",chi_dof)
    print("       pcov: \n",pcov,"\n")

    fitter.plot_fit()
    plt.show()


if __name__ == '__main__':
    fitExample()