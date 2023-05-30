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


def continuumExtrapolate(x,obs,obs_err,order=1,show_results=False,plot_results=False,paramLabels=None,prior=None,
                         prior_err=None,error_strat='propagation',algorithms=None,xtype="a",plotName="contExtrap.pdf",**kwargs):
    """ Do a continuum limit extrapolation at some order in a^2. Allows the option for priors in case you want
        to fit to a higher order series, and you have some idea what the coefficients should be like. """
    if order<1:
        logger.TBError('Please input order > 1.')

    if algorithms is None:
        algorithms = ["L-BFGS-B", "TNC", "Powell", "Nelder-Mead", "COBYLA", "CG", "BFGS", "dogleg", "trust-ncg"]

    coeffs = ()
    for i in range(order+1):
        coeffs += (1.0,)

    if xtype == "a":
        x = np.array(x)**2
    elif xtype == "Nt":
        x = 1/np.array(x)**2
    else:
        logger.TBError('Unknown xtype',xtype)

    fit = Fitter(powerSeries, x, obs, obs_err, norm_err_chi2=False, func_sup_numpy=False, expand=False,
                 error_strat = error_strat)

    if prior is None:
        result, result_err, chidof = fit.try_fit(start_params=coeffs, priorval=prior, priorsigma=prior_err, algorithms=algorithms)
    else:
        result, result_err, chidof, logGBF, _ = fit.try_fit(start_params=coeffs, priorval=prior, priorsigma=prior_err,
                                                            algorithms=algorithms, detailedInfo=True)

    if show_results:
        logger.info()
        for i in range(len(coeffs)):
            if paramLabels is None:
                logger.info('        c_'+str(i)+' = '+get_err_str(result[i],result_err[i]))
            else:
                logger.info(paramLabels[i] + ' = ' + get_err_str(result[i], result_err[i]))
        logger.info('chi2/d.o.f. =',round(chidof,3))
        if prior is not None:
            logger.info('     logGBF =', round(logGBF, 3))
        logger.info()

    if plot_results:
        latexify()
        fit.plot_fit(**kwargs)
        plot_dots([0],result[0],result_err[0],color='red')
        plt.savefig(plotName)
        plt.show()

    return result, result_err, chidof