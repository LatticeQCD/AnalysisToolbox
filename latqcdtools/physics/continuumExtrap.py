# 
# continuumExtrap.py                                                               
# 
# D. Clarke
# 
# A simple module for performing extrapolations to fit functions.
# 

import numpy as np
from latqcdtools.base.printErrorBars import get_err_str
import latqcdtools.base.logger as logger
from latqcdtools.base.plotting import plt, plot_dots, latexify
from latqcdtools.statistics.fitting import Fitter


def powerSeries(x,coeffs):
    result = 0.
    for i in range(len(coeffs)):
        result += coeffs[i]*x**i
    return result


def extrapolate_from_a(a,obs,obs_err,order=1,show_results=False,plot_results=False,obsName=None):
    """ Do a continuum limit extrapolation at some order in a^2. """
    if order<1:
        logger.TBError('Please input order > 1.')
    coeffs = ()
    for i in range(order+1):
        coeffs += (1.0,)
    a = np.array(a)**2
    fit = Fitter( powerSeries, a, obs, obs_err, func_sup_numpy=False, expand=False )
    result, result_err, chidof = fit.do_fit(algorithm = 'curve_fit',start_params=coeffs)
    if show_results:
        for i in range(len(coeffs)):
            print('c_'+str(i)+' = '+get_err_str(result[i],result_err[i]))
        print('chi2/d.o.f. = ',round(chidof,3))
    if plot_results:
        latexify()
        fit.plot_fit(xlabel='$a^2$ [fm]',ylabel=obsName)
        plot_dots([0],result[0],result_err[0],color='red')
        plt.show()
    return result, result_err, chidof