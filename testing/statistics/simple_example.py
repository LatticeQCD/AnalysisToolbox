# 
# simple_example.py                                                               
# 
# H. Sandmeyer
# 
# Example of a simple 3-parameter quadratic fit.
#
import numpy as np
import matplotlib.pyplot as plt
from latqcdtools.statistics.fitting import Fitter

print("\n Example of a simple 3-parameter quadratic fit.\n")

def fit_func(x, params):
    a = params[0]
    b = params[1]
    c = params[2]
    return a*x**2 + b*x + c

xdata, ydata, edata = np.genfromtxt("wurf.dat", usecols=(0,2,3), unpack=True)

fitter = Fitter(fit_func, xdata, ydata, expand = False)

res, res_err, chi_dof, pcov = fitter.try_fit(start_params = [1, 2, 3], algorithms = ['levenberg', 'curve_fit'], 
                                             ret_pcov = True)

print(" a , b,  c : ",res)
print(" ae, be, ce: ",res_err)
print("chi2/d.o.f.: ",chi_dof)
print("       pcov: \n",pcov,"\n")

fitter.plot_fit()
plt.show()
