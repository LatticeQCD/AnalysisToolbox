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
from latqcdtools.base.plotting import plt, plot_dots, fill_param_dict
from latqcdtools.statistics.fitting import Fitter, std_algs, bayes_algs
from latqcdtools.base.speedify import DEFAULTTHREADS
from latqcdtools.base.check import checkType


def powerSeries(x,coeffs):
    result = 0.
    for i in range(len(coeffs)):
        result += coeffs[i]*x**i
    return result


class Extrapolator(Fitter):

    def __init__(self, x, obs, obs_err, order=1, xtype="a", error_strat='propagation', nproc=DEFAULTTHREADS):
        """ A framework for doing continuum limit extrapolations. Assume a power series to some order in a^2.

        Args:
            x (array-like): a data or Nt data 
            obs (array-like)
            obs_err (array-like)
            order (int, optional): order of the power series. Defaults to 1.
            xtype (str, optional): choose to fit a data or Nt data. Defaults to "a".
            error_strat (str, optional): calculate errors using error propagation or augmented chi^2. Defaults to 'propagation'.
            nproc (int, optional): number of processors for fitting. Defaults to DEFAULTTHREADS.
        """
        checkType(order,int)

        self._order              = order
        self._triedExtrapolation = False
        self._logGBF             = None

        if xtype == "a":
            x = np.array(x)**2
        elif xtype == "Nt":
            x = 1/np.array(x)**2
        else:
            logger.TBError('Unknown xtype',xtype)
        if order<1:
            logger.TBError('Please input order > 1.')

        Fitter.__init__(self, powerSeries, x, obs, obs_err, norm_err_chi2=False, expand=False, error_strat=error_strat, nproc=nproc)

    def extrapolate(self,start_coeffs=None,prior=None,prior_err=None):
        """ Carry out the extrapolation.

        Args:
            start_coeffs (array-like, optional): your guess for starting parameters. Defaults to None.
            prior (array-like, optional): Bayesian priors. Defaults to None.
            prior_err (array-like, optional): Bayesian prior errors. Defaults to None.

        Returns:
            (array-like, array-like, float, float): fit result, its error, chi^2/d.o.f., and logGBF if relevant.
        """
        if start_coeffs is None:
            coeffs = ()
            for i in range(self._order+1):
                coeffs += (1.0,)
        else:
            coeffs=start_coeffs
        self._triedExtrapolation = True
        if prior is None:
            self._result, self._result_err, self._chidof = self.try_fit(start_params=coeffs, algorithms=std_algs)
            return self._result, self._result_err, self._chidof
        else:
            self._result, self._result_err, self._chidof, self._logGBF, _ = self.try_fit(start_params=coeffs, priorval=prior, priorsigma=prior_err, 
                                                                                         algorithms=bayes_algs, detailedInfo=True)
            return self._result, self._result_err, self._chidof, self._logGBF

    def plot(self,**kwargs):
        """ Add extrapolation to plot. Accepts the same kwargs as Fitter.plot_fit. """
        if not self._triedExtrapolation:
            logger.TBError("Can't plot an extrapolation without having extrapolated first...")
        plot_dots([0],self._result[0],self._result_err[0],color=kwargs['color'])
        self.plot_data(color=kwargs['color'])
        self.plot_fit(**kwargs)

    def showResults(self):
        """ Print extrapolation results to screen. """
        if not self._triedExtrapolation:
            logger.TBError("Can't show extrapolation results without having extrapolated first...")
        logger.info()
        for i in range(len(self._result)):
            logger.info('        c_'+str(i)+' = '+get_err_str(self._result[i],self._result_err[i]))
        logger.info('chi2/d.o.f. =',round(self._chidof,3))
        if self._logGBF is not None:
            logger.info('     logGBF =',round(self._logGBF, 3))
        logger.info()


def continuumExtrapolate(x,obs,obs_err,order=1,show_results=False,plot_results=False,prior=None, start_coeffs=None,prior_err=None,
                         error_strat='propagation',xtype="a",nproc=DEFAULTTHREADS):
    """ A convenience wrapper for the Extrapolator. """
    ext = Extrapolator(x, obs, obs_err, order=order, error_strat=error_strat, nproc=nproc)
    result = ext.extrapolate(start_coeffs=start_coeffs, prior=prior, prior_err=prior_err)
    if show_results:
        ext.showResults()
    if plot_results:
        ext.plot()
        plt.show()
    return result