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
from latqcdtools.base.plotting import plt
from latqcdtools.statistics.fitting import Fitter, std_algs, bayes_algs
from latqcdtools.base.speedify import DEFAULTTHREADS
from latqcdtools.base.check import checkType


def _powerSeries(x,coeffs):
    """ The default fit form for a continuum extrapolation is a power series in a^2.

    Args:
        x (array-like): a^2 or 1/Nt^2 
        coeffs (array-like): power series coefficients 

    Returns:
        array-like: power series in x 
    """
    result = 0.
    for i in range(len(coeffs)):
        result += coeffs[i]*x**i
    return result


class Extrapolator(Fitter):


    def __init__(self, x, obs, obs_err, order=1, xtype="a", error_strat='propagation', ansatz=None, nproc=DEFAULTTHREADS):
        """ A framework for doing continuum limit extrapolations.

        Args:
            x (array-like): a^2 data or 1/Nt^2 data 
            obs (array-like)
            obs_err (array-like)
            ansatz (func, optional): continuum-limit fit ansatz. Power series in a^2 by default.
            order (int, optional): order of the power series. Defaults to 1.
            xtype (str, optional): choose to fit a data or Nt data. Defaults to "a".
            error_strat (str, optional): calculate errors using error propagation or augmented chi^2. Defaults to 'propagation'.
            nproc (int, optional): number of processors for fitting. Defaults to DEFAULTTHREADS.
        """
        checkType(order,int)

        self._order              = order
        self._triedExtrapolation = False
        self._ansatz             = ansatz

        if xtype == "a":
            x = np.array(x)**2
        elif xtype == "Nt":
            x = 1/np.array(x)**2
        else:
            logger.TBError('Unknown xtype',xtype)
        if order<1:
            logger.TBError('Please input order > 1.')

        if ansatz is None:
            ansatz = _powerSeries
        else:
            if self._order != 1:
                logger.warn('Not using a power series ansatz, but still using custom order.')
            logger.debug('received ansatz',ansatz)

        Fitter.__init__(self, ansatz, x, obs, obs_err, norm_err_chi2=False, error_strat=error_strat, nproc=nproc)


    def __repr__(self) -> str:
        return "Extrapolator"


    def extrapolate(self,start_coeffs=None,prior=None,priorsigma=None,detailedInfo=False):
        """ Carry out the extrapolation.

        Args:
            start_coeffs (array-like, optional): your guess for starting parameters. Defaults to None.
            prior (array-like, optional): Bayesian priors. Defaults to None.
            priorsigma (array-like, optional): Bayesian prior errors. Defaults to None.
            detailedInfo (bool, optional): Do you want information like AIC? Defaults to False.
            
        Returns:
            (array-like, array-like, float, float): fit result, its error, chi^2/d.o.f., and logGBF if relevant.
        """
        if start_coeffs is None:
            if self._ansatz is not None:
                logger.TBError('You need to provide start_coeffs if you use an ansatz.')
            coeffs = ()
            for i in range(self._order+1):
                coeffs += (1.0,)
        else:
            coeffs=start_coeffs
        self._triedExtrapolation = True
        if prior is None:
            self._result, self._result_err, self._chidof, self._stats = self.try_fit(start_params=coeffs, algorithms=std_algs, 
                                                                                     detailedInfo=True)
        else:
            self._result, self._result_err, self._chidof, self._stats = self.try_fit(start_params=coeffs, priorval=prior, 
                                                                                     priorsigma=priorsigma, 
                                                                                     algorithms=bayes_algs, detailedInfo=True)
        if detailedInfo:
            return self._result, self._result_err, self._chidof, self._stats
        else:
            return self._result, self._result_err, self._chidof


    def plot(self,**kwargs):
        """ Add extrapolation to plot. Accepts the same kwargs as Fitter.plot_fit. """
        if not self._triedExtrapolation:
            logger.TBError("Can't plot an extrapolation without having extrapolated first...")
        domain=(1e-8,np.max(self._xdata))
        self.plot_data(**kwargs)
        kwargs['label']=None
        self.plot_fit(domain,**kwargs)


    def save(self,filename,header):
        if not self._triedExtrapolation:
            logger.TBError("Can't save an extrapolation without having extrapolated first...")
        domain=(1e-8,np.max(self._xdata))
        self.save_func(filename,domain,header=header)


    def showResults(self):
        """ Print extrapolation results to screen. """
        if not self._triedExtrapolation:
            logger.TBError("Can't show extrapolation results without having extrapolated first...")
        logger.info()
        for i in range(len(self._result)):
            logger.info('        c_'+str(i)+' = '+get_err_str(self._result[i],self._result_err[i]))
        logger.info('chi2/d.o.f. =',round(self._chidof,3))
        logger.info('     logGBF =',round(self._stats['logGBF'], 3))
        logger.info()


def continuumExtrapolate(x,obs,obs_err,order=1,show_results=False,plot_results=False,prior=None, start_coeffs=None,priorsigma=None,
                         error_strat='propagation',xtype="a",nproc=DEFAULTTHREADS,detailedInfo=False):
    """ A convenience wrapper for the Extrapolator. """
    ext = Extrapolator(x, obs, obs_err, xtype=xtype, order=order, error_strat=error_strat, nproc=nproc)
    result = ext.extrapolate(start_coeffs=start_coeffs, prior=prior, priorsigma=priorsigma,detailedInfo=detailedInfo)
    if show_results:
        ext.showResults()
    if plot_results:
        if xtype=="a":
            xlabel = "$a^2$"
        else:
            xlabel = "$1/N_\\tau^2$"
        ext.plot(xlabel=xlabel)
        plt.show()
    return result
