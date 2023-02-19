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


def extrapolate_from_a(a,obs,obs_err,order=1,show_results=False,plot_results=False,paramLabels=None,prior=None,
                       prior_err=None,**kwargs):
    """ Do a continuum limit extrapolation at some order in a^2. Allows the option for priors in case you want
        to fit to a higher order series, and you have some idea what the coefficients should be like. """
    if order<1:
        logger.TBError('Please input order > 1.')

    coeffs = ()
    for i in range(order+1):
        coeffs += (1.0,)
    a = np.array(a)**2
    fit = Fitter( powerSeries, a, obs, obs_err, norm_err_chi2=False, cut_eig=True, func_sup_numpy=False, expand=False )
    result, result_err, chidof = fit.try_fit(start_params=coeffs, priorval=prior, priorsigma=prior_err)

    if show_results:
        print()
        for i in range(len(coeffs)):
            if paramLabels is None:
                print('c_'+str(i)+' = '+get_err_str(result[i],result_err[i]))
            else:
                print(paramLabels[i] + ' = ' + get_err_str(result[i], result_err[i]))
        print('chi2/d.o.f. = ',round(chidof,3),'\n')

    if plot_results:
        latexify()
        fit.plot_fit(**kwargs)
        plot_dots([0],result[0],result_err[0],color='red')
        plt.show()

    return result, result_err, chidof